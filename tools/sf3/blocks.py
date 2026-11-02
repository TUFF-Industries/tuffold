import sys, struct
path, size, start, count = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
data = open(path,'rb').read()
for b in range(start, start+count):
    blk = data[b*size:(b+1)*size]
    if len(blk) < size: break
    words = [blk[i:i+4] for i in range(0, size - size % 4, 4)]
    ints = [struct.unpack('<i', w)[0] for w in words]
    fl = [struct.unpack('<f', w)[0] for w in words]
    cells = []
    for i,(n,f) in enumerate(zip(ints, fl)):
        if -100000 < n < 1000000: cells.append(f"{i*4}:{n}")
        elif abs(f) > 1e-6 and abs(f) < 1e7: cells.append(f"{i*4}:{f:.5g}f")
        else: cells.append(f"{i*4}:{words[i].hex()}")
    print(f"#{b}", " ".join(cells), ("tail:"+blk[size - size%4:].hex()) if size%4 else "")
