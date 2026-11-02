"""Reader for Scaffold's block-based object records (hypothesis, checked by
parsing every record to its end).

A record: int32 LE id, int32 LE (−1 = no continuation?), then objects in a
Java-serialization dialect: 0x73 object / 0x75 array / 0x70 null / 0x74
string, class references as 0x71 + int32 LE class id (1-based index into the
descriptor's class list), array lengths int32 LE, field values big-endian in
class-descriptor order (superclass first). Classes with writeObject (flag
0x01) are followed by custom data until 0x78.
"""
import struct
from jser import ClassDesc


class BBObj:
    def __init__(self, cls):
        self.cls, self.fields, self.custom = cls, {}, []

    def __repr__(self):
        return f"BBObj({self.cls.name.rsplit('.', 1)[-1]})"


class BBReader:
    def __init__(self, data, classes, pos=0):
        self.d, self.p, self.classes = data, pos, classes

    def read(self, n):
        b = self.d[self.p:self.p + n]
        if len(b) != n:
            raise EOFError(self.p)
        self.p += n
        return b

    def cls(self):
        tc = self.read(1)[0]
        if tc == 0x70:
            return None
        if tc != 0x71:
            raise ValueError(f"class ref expected at {self.p - 1}: {tc:#x}")
        cid = struct.unpack("<i", self.read(4))[0]
        return self.classes[cid - 1]

    def value(self, code):
        if code == "B": return struct.unpack(">b", self.read(1))[0]
        if code == "C": return chr(struct.unpack(">H", self.read(2))[0])
        if code == "D": return struct.unpack(">d", self.read(8))[0]
        if code == "F": return struct.unpack(">f", self.read(4))[0]
        if code == "I": return struct.unpack(">i", self.read(4))[0]
        if code == "J": return struct.unpack(">q", self.read(8))[0]
        if code == "S": return struct.unpack(">h", self.read(2))[0]
        if code == "Z": return self.read(1)[0] != 0
        return self.obj()

    def obj(self):
        tc = self.read(1)[0]
        if tc == 0x70:
            return None
        if tc == 0x74:
            n = struct.unpack("<H", self.read(2))[0]  # LE, unlike Java's
            return self.read(n).decode("utf-8", "replace")
        if tc == 0x75:
            c = self.cls()
            n = struct.unpack("<i", self.read(4))[0]
            code = c.name[1]
            sizes = {"B": 1, "Z": 1, "S": 2, "C": 2, "I": 4, "F": 4, "J": 8, "D": 8}
            if code in sizes:
                fmt = {"B": "b", "Z": "?", "S": "h", "C": "H", "I": "i", "F": "f", "J": "q", "D": "d"}[code]
                return list(struct.unpack(f"<{n}{fmt}", self.read(sizes[code] * n)))
            return [self.obj() for _ in range(n)]
        if tc == 0x73:
            c = self.cls()
            o = BBObj(c)
            for desc in c.chain():
                for code, name, _ in desc.fields:
                    o.fields[name] = self.value(code)
                if desc.flags & 0x01:
                    o.custom.append(self.custom(desc.name))
            return o
        if tc == 0x7E:  # enum constant: class, then its name as a string
            self.cls()
            return self.obj()
        raise ValueError(f"unknown tag {tc:#x} at {self.p - 1}")

    def i4le(self):
        return struct.unpack("<i", self.read(4))[0]

    def custom(self, name):
        """writeObject data, little-endian, no end marker: per class."""
        if name == "java.util.TreeMap" or name == "java.util.HashMap" and False:
            n = self.i4le()
            return [(self.obj(), self.obj()) for _ in range(n)]
        if name == "java.util.HashSet":  # capacity, load factor, size, elements
            self.i4le(); self.read(4); n = self.i4le()
            return [self.obj() for _ in range(n)]
        if name == "gnu.trove.TIntArrayList":
            n = self.i4le()
            return [self.i4le() for _ in range(n)]
        raise ValueError(f"no custom reader for {name} at {self.p}")


def class_list(descriptor):
    cm = descriptor.fields["serialization"].fields["cm"]
    return [x for a in cm.fields["classPos"].annotations.values() for x in a if isinstance(x, ClassDesc)]
