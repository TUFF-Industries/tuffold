import sys, struct
def read_strings(path, block=68):
    data = open(path, 'rb').read()
    n = len(data) // block
    blocks = [data[i*block:(i+1)*block] for i in range(n)]
    heads = {}
    cont = set()
    for i, b in enumerate(blocks):
        nxt = struct.unpack('<i', b[-4:])[0]
        if nxt not in (-1,) and 0 < nxt < n: cont.add(nxt)
    out = {}
    for i, b in enumerate(blocks):
        if i == 0 or i in cont: continue
        rid = struct.unpack('<i', b[:4])[0]
        if rid <= 0: continue
        text, j, seen = b'', i, 0
        while True:
            blk = blocks[j]
            payload = blk[4:-4] if j == i else blk[4:-4]
            text += payload
            nxt = struct.unpack('<i', blk[-4:])[0]
            seen += 1
            if nxt == -1 or nxt <= 0 or nxt >= n or seen > 1000: break
            j = nxt
        out[rid] = text.split(b'\x02')[0].decode('latin-1')
    return out
if __name__ == '__main__':
    s = read_strings(sys.argv[1])
    print(len(s), 'records; ids', min(s), '..', max(s))
    for k in list(s)[:5] + list(s)[-3:]: print(k, repr(s[k][:120]))
