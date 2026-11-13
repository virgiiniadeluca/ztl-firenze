import struct, json, math, sys, os

# --- DBF ---
def read_dbf(p):
    b = open(p, 'rb').read()
    n = struct.unpack('<I', b[4:8])[0]
    hl, rl = struct.unpack('<HH', b[8:12])
    fields = []
    o = 32
    while b[o] != 0x0d:
        fields.append((b[o:o+11].split(b'\0')[0].decode(), chr(b[o+11]), b[o+16]))
        o += 32
    rows = []
    for i in range(n):
        r = b[hl+i*rl:hl+(i+1)*rl]
        pos = 1
        d = {}
        for name, t, l in fields:
            raw = r[pos:pos+l]
            pos += l
            try:
                s = raw.decode('utf-8')
            except UnicodeDecodeError:
                s = raw.decode('latin-1')
            s = s.strip()
            if t in 'NF' and s:
                try:
                    s = float(s) if ('.' in s) else int(s)
                except ValueError:
                    pass
            d[name] = s
        rows.append(d)
    return rows

# --- SHP ---
def read_shp(p):
    b = open(p, 'rb').read()
    stype = struct.unpack('<i', b[32:36])[0]
    o = 100
    out = []
    while o < len(b):
        clen = struct.unpack('>i', b[o+4:o+8])[0] * 2
        rec = b[o+8:o+8+clen]
        o += 8 + clen
        t = struct.unpack('<i', rec[:4])[0]
        if t == 0:
            out.append(None); continue
        if t == 1:
            x, y = struct.unpack('<2d', rec[4:20]); out.append((x, y))
        elif t == 5:
            nparts, npts = struct.unpack('<2i', rec[36:44])
            parts = list(struct.unpack('<%di' % nparts, rec[44:44+4*nparts]))
            q = 44 + 4*nparts
            pts = [struct.unpack('<2d', rec[q+16*i:q+16*i+16]) for i in range(npts)]
            parts.append(npts)
            out.append([pts[parts[i]:parts[i+1]] for i in range(nparts)])
        else:
            raise SystemExit('tipo non gestito %d' % t)
    return stype, out

# --- Monte Mario / Italy zone 1 (EPSG:3003) -> WGS84 ---
a_i, f_i = 6378388.0, 1/297.0
a_w, f_w = 6378137.0, 1/298.257223563
TX, TY, TZ = -104.1, -49.1, -9.9
RX, RY, RZ = [math.radians(v/3600.0) for v in (0.971, -2.917, 0.714)]
DS = -11.68e-6
LON0 = math.radians(9.0); K0 = 0.9996; FE = 1500000.0; FN = 0.0

def tm_inverse(E, N, a, f):
    e2 = f*(2-f); ep2 = e2/(1-e2)
    x = E - FE; y = N - FN
    M = y / K0
    mu = M / (a*(1 - e2/4 - 3*e2**2/64 - 5*e2**3/256))
    e1 = (1 - math.sqrt(1-e2)) / (1 + math.sqrt(1-e2))
    phi1 = (mu + (3*e1/2 - 27*e1**3/32)*math.sin(2*mu)
            + (21*e1**2/16 - 55*e1**4/32)*math.sin(4*mu)
            + (151*e1**3/96)*math.sin(6*mu)
            + (1097*e1**4/512)*math.sin(8*mu))
    s, c, t = math.sin(phi1), math.cos(phi1), math.tan(phi1)
    N1 = a / math.sqrt(1 - e2*s*s)
    T1 = t*t; C1 = ep2*c*c
    R1 = a*(1-e2) / (1 - e2*s*s)**1.5
    D = x / (N1*K0)
    phi = phi1 - (N1*t/R1)*(D**2/2 - (5 + 3*T1 + 10*C1 - 4*C1**2 - 9*ep2)*D**4/24
          + (61 + 90*T1 + 298*C1 + 45*T1**2 - 252*ep2 - 3*C1**2)*D**6/720)
    lam = LON0 + (D - (1 + 2*T1 + C1)*D**3/6
          + (5 - 2*C1 + 28*T1 - 3*C1**2 + 8*ep2 + 24*T1**2)*D**5/120) / c
    return phi, lam

def geo2xyz(phi, lam, a, f):
    e2 = f*(2-f); N = a/math.sqrt(1-e2*math.sin(phi)**2)
    return (N*math.cos(phi)*math.cos(lam), N*math.cos(phi)*math.sin(lam), N*(1-e2)*math.sin(phi))

def xyz2geo(X, Y, Z, a, f):
    e2 = f*(2-f); lam = math.atan2(Y, X); p = math.hypot(X, Y)
    phi = math.atan2(Z, p*(1-e2))
    for _ in range(10):
        N = a/math.sqrt(1-e2*math.sin(phi)**2)
        phi = math.atan2(Z + e2*N*math.sin(phi), p)
    return phi, lam

def gb_to_wgs84(E, N):
    phi, lam = tm_inverse(E, N, a_i, f_i)
    X, Y, Z = geo2xyz(phi, lam, a_i, f_i)
    # position vector rotation (EPSG 9606)
    X2 = TX + (1+DS)*(X - RZ*Y + RY*Z)
    Y2 = TY + (1+DS)*(RZ*X + Y - RX*Z)
    Z2 = TZ + (1+DS)*(-RY*X + RX*Y + Z)
    phi2, lam2 = xyz2geo(X2, Y2, Z2, a_w, f_w)
    return math.degrees(lam2), math.degrees(phi2)

def conv(pt):
    lon, lat = gb_to_wgs84(*pt)
    return [round(lon, 7), round(lat, 7)]

def area(r):
    return sum(r[i][0]*r[i+1][1] - r[i+1][0]*r[i][1] for i in range(len(r)-1))/2

def pip(pt, ring):
    x, y = pt; ins = False
    for i in range(len(ring)-1):
        x1, y1 = ring[i]; x2, y2 = ring[i+1]
        if (y1 > y) != (y2 > y) and x < (x2-x1)*(y-y1)/(y2-y1)+x1:
            ins = not ins
    return ins

def polygon_geom(rings):
    outers, holes = [], []
    for r in rings:
        (outers if area(r) < 0 else holes).append(r)  # shapefile: outer = clockwise
    polys = [[o] for o in outers]
    for h in holes:
        for pl in polys:
            if pip(h[0], pl[0]):
                pl.append(h); break
    res = []
    for pl in polys:
        rr = []
        for i, r in enumerate(pl):
            c = [conv(p) for p in r]
            # GeoJSON RFC 7946: outer CCW, holes CW
            wantccw = (i == 0)
            if (area(c) > 0) != wantccw:
                c.reverse()
            rr.append(c)
        res.append(rr)
    return {'type': 'Polygon', 'coordinates': res[0]} if len(res) == 1 else {'type': 'MultiPolygon', 'coordinates': res}

def convert(base, out, validate=None):
    stype, geoms = read_shp(base + '.shp')
    rows = read_dbf(base + '.dbf')
    assert len(geoms) == len(rows), (len(geoms), len(rows))
    feats = []
    for g, r in zip(geoms, rows):
        if g is None: continue
        geom = ({'type': 'Point', 'coordinates': conv(g)} if stype == 1 else polygon_geom(g))
        feats.append({'type': 'Feature', 'properties': r, 'geometry': geom})
    json.dump({'type': 'FeatureCollection', 'features': feats}, open(out, 'w', encoding='utf-8'), ensure_ascii=False, separators=(',', ':'))
    print(out, len(feats), 'elementi', os.path.getsize(out), 'byte')
    return feats

if __name__ == '__main__':
    src, dst = sys.argv[1], sys.argv[2]
    os.makedirs(dst, exist_ok=True)
    f = convert(src + '/varchi_scudoverde/varchi_scudoverde', dst + '/varchi_scudo_verde.geojson')
    # controllo di precisione: coordinate GPS ufficiali del dataset vs conversione
    errs = []
    for ft in f:
        p = ft['properties']; lon, lat = ft['geometry']['coordinates']
        dx = (lon - float(p['coord_x'])) * 111320 * math.cos(math.radians(lat))
        dy = (lat - float(p['coord_y'])) * 110574
        errs.append(math.hypot(dx, dy))
    print('Errore conversione vs coordinate ufficiali: medio %.2f m, massimo %.2f m su %d varchi' % (sum(errs)/len(errs), max(errs), len(errs)))
    convert(src + '/ztl_diurna/ztl_diurnaPolygon', dst + '/ztl_diurna.geojson')
    convert(src + '/ztl_notturna/ztl_notturnaPolygon', dst + '/ztl_notturna.geojson')
    convert(src + '/sas_porte_telematiche_od/sas_porte_telematiche_od', dst + '/porte_telematiche.geojson')
