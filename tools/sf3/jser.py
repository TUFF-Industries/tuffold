"""Minimal Java Object Serialization Stream reader (public spec, protocol 2).

Parses into plain Python structures:
  Obj(cls, fields: dict, annotations: list)   # annotations = block data / objects written by writeObject
  Arr(cls, items)
  str, None, Ref handled by sharing the same Python object.
Never executes anything; it only reads bytes.
"""
import struct

TC_NULL, TC_REFERENCE, TC_CLASSDESC, TC_OBJECT, TC_STRING, TC_ARRAY, TC_CLASS = 0x70, 0x71, 0x72, 0x73, 0x74, 0x75, 0x76
TC_BLOCKDATA, TC_ENDBLOCKDATA, TC_RESET, TC_BLOCKDATALONG, TC_EXCEPTION = 0x77, 0x78, 0x79, 0x7A, 0x7B
TC_LONGSTRING, TC_PROXYCLASSDESC, TC_ENUM = 0x7C, 0x7D, 0x7E
BASE = 0x7E0000
SC_WRITE_METHOD, SC_SERIALIZABLE, SC_EXTERNALIZABLE, SC_BLOCK_DATA, SC_ENUM = 0x01, 0x02, 0x04, 0x08, 0x10


class ClassDesc:
    def __init__(self, name, suid, flags, fields, superdesc):
        self.name, self.suid, self.flags, self.fields, self.superdesc = name, suid, flags, fields, superdesc

    def chain(self):
        out, c = [], self
        while c is not None:
            out.append(c)
            c = c.superdesc
        return list(reversed(out))

    def __repr__(self):
        return f"ClassDesc({self.name})"


class Obj:
    def __init__(self, cls):
        self.cls, self.fields, self.annotations = cls, {}, {}

    def __repr__(self):
        return f"Obj({self.cls.name})"


class Arr:
    def __init__(self, cls, items):
        self.cls, self.items = cls, items

    def __repr__(self):
        return f"Arr({self.cls.name}, {len(self.items)})"


class Block(bytes):
    pass


class Enum:
    def __init__(self, cls, name):
        self.cls, self.name = cls, name

    def __repr__(self):
        return f"Enum({self.cls.name}.{self.name})"


class Reader:
    def __init__(self, data):
        self.d, self.p, self.handles = data, 0, []

    def read(self, n):
        b = self.d[self.p:self.p + n]
        if len(b) != n:
            raise EOFError(f"wanted {n} at {self.p}")
        self.p += n
        return b

    def u1(self): return self.read(1)[0]
    def u2(self): return struct.unpack(">H", self.read(2))[0]
    def i4(self): return struct.unpack(">i", self.read(4))[0]
    def i8(self): return struct.unpack(">q", self.read(8))[0]

    def utf(self):
        return self.read(self.u2()).decode("utf-8", "replace")

    def new_handle(self, value):
        self.handles.append(value)
        return len(self.handles) - 1

    def stream(self):
        magic, version = self.u2(), self.u2()
        assert magic == 0xACED and version == 5, (hex(magic), version)
        out = []
        while self.p < len(self.d):
            out.append(self.content())
        return out

    def content(self):
        tc = self.u1()
        if tc == TC_BLOCKDATA:
            return Block(self.read(self.u1()))
        if tc == TC_BLOCKDATALONG:
            return Block(self.read(self.i4()))
        return self.object(tc)

    def object(self, tc=None):
        if tc is None:
            tc = self.u1()
        if tc == TC_NULL:
            return None
        if tc == TC_REFERENCE:
            return self.handles[self.i4() - BASE]
        if tc == TC_STRING:
            s = self.utf()
            self.new_handle(s)
            return s
        if tc == TC_LONGSTRING:
            s = self.read(self.i8()).decode("utf-8", "replace")
            self.new_handle(s)
            return s
        if tc in (TC_CLASSDESC, TC_PROXYCLASSDESC):
            return self.classdesc(tc)
        if tc == TC_CLASS:
            value = ("class", self.classdesc())
            self.new_handle(value)
            return value
        if tc == TC_OBJECT:
            return self.new_object()
        if tc == TC_ARRAY:
            return self.array()
        if tc == TC_ENUM:
            cls = self.classdesc()
            e = Enum(cls, None)
            self.new_handle(e)
            e.name = self.object()
            return e
        if tc == TC_RESET:
            self.handles = []
            return self.object()
        raise ValueError(f"unknown tc {tc:#x} at {self.p - 1}")

    def classdesc(self, tc=None):
        if tc is None:
            tc = self.u1()
        if tc == TC_NULL:
            return None
        if tc == TC_REFERENCE:
            return self.handles[self.i4() - BASE]
        if tc == TC_PROXYCLASSDESC:
            n = self.i4()
            names = [self.utf() for _ in range(n)]
            desc = ClassDesc("Proxy" + str(names), 0, SC_SERIALIZABLE, [], None)
            self.new_handle(desc)
            self.annotation()
            desc.superdesc = self.classdesc()
            return desc
        assert tc == TC_CLASSDESC, hex(tc)
        name = self.utf()
        suid = self.i8()
        desc = ClassDesc(name, suid, 0, [], None)
        self.new_handle(desc)
        desc.flags = self.u1()
        for _ in range(self.u2()):
            code = chr(self.u1())
            fname = self.utf()
            ftype = self.object() if code in "[L" else code
            desc.fields.append((code, fname, ftype))
        self.annotation()
        desc.superdesc = self.classdesc()
        return desc

    def annotation(self):
        out = []
        while True:
            if self.d[self.p] == TC_ENDBLOCKDATA:
                self.p += 1
                return out
            out.append(self.content())

    def value(self, code):
        if code == "B": return struct.unpack(">b", self.read(1))[0]
        if code == "C": return chr(self.u2())
        if code == "D": return struct.unpack(">d", self.read(8))[0]
        if code == "F": return struct.unpack(">f", self.read(4))[0]
        if code == "I": return self.i4()
        if code == "J": return self.i8()
        if code == "S": return struct.unpack(">h", self.read(2))[0]
        if code == "Z": return self.u1() != 0
        return self.object()

    def new_object(self):
        cls = self.classdesc()
        obj = Obj(cls)
        self.new_handle(obj)
        for desc in cls.chain():
            if desc.flags & SC_EXTERNALIZABLE:
                obj.annotations[desc.name] = self.annotation()
                continue
            for code, fname, _ in desc.fields:
                obj.fields[fname] = self.value(code)
            if desc.flags & SC_WRITE_METHOD:
                obj.annotations[desc.name] = self.annotation()
        return obj

    def array(self):
        cls = self.classdesc()
        n = self.i4()
        arr = Arr(cls, [])
        self.new_handle(arr)
        code = cls.name[1]
        if code in "BZ":
            arr.items = bytes(self.read(n))
        elif code == "I":
            arr.items = list(struct.unpack(f">{n}i", self.read(4 * n)))
        elif code == "J":
            arr.items = list(struct.unpack(f">{n}q", self.read(8 * n)))
        elif code == "D":
            arr.items = list(struct.unpack(f">{n}d", self.read(8 * n)))
        elif code == "F":
            arr.items = list(struct.unpack(f">{n}f", self.read(4 * n)))
        elif code == "S":
            arr.items = list(struct.unpack(f">{n}h", self.read(2 * n)))
        elif code == "C":
            arr.items = "".join(chr(c) for c in struct.unpack(f">{n}H", self.read(2 * n)))
        else:
            arr.items = [self.object() for _ in range(n)]
        return arr


def load(path):
    return Reader(open(path, "rb").read()).stream()


def short(name):
    return name.rsplit(".", 1)[-1] if name else name
