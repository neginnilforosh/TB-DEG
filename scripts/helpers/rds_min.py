"""Minimal reader for R .rds files (XDR, format version 2/3). Read-only; supports the SEXP types found in
WGCNA workspace objects (vectors, lists, pairlists, symbols, language objects, attributes, references).
Written for the TB-DEG pre-TB-HIV checkpoint because no R installation is available."""
import gzip, struct, bz2, lzma

class RObj:
    def __init__(self, t, value=None, attrs=None):
        self.type, self.value, self.attrs = t, value, attrs or {}
    def attr(self, k, default=None):
        return self.attrs.get(k, default)
    def __repr__(self):
        return f"RObj({self.type}, attrs={list(self.attrs)})"

NA_INT = -2147483648

class Reader:
    def __init__(self, data):
        self.d, self.p, self.refs = data, 0, []
    def i(self):
        v = struct.unpack(">i", self.d[self.p:self.p + 4])[0]; self.p += 4; return v
    def raw(self, n):
        b = self.d[self.p:self.p + n]; self.p += n; return b
    def header(self):
        assert self.raw(2) == b"X\n", "only XDR format supported"
        ver = self.i(); self.i(); self.i()
        if ver == 3:
            n = self.i(); self.raw(n)
    def item(self):
        flags = self.i()
        t = flags & 0xFF; has_attr = bool(flags & (1 << 9)); has_tag = bool(flags & (1 << 10))
        if t == 254: return None                        # NILVALUE_SXP
        if t == 242: return RObj("emptyenv")
        if t == 241: return RObj("baseenv")
        if t == 253: return RObj("globalenv")
        if t == 252: return RObj("unboundvalue")
        if t == 251: return RObj("missingarg")
        if t == 250: return RObj("basenamespace")
        if t == 255:                                    # REFSXP
            idx = flags >> 8
            if idx == 0: idx = self.i()
            return self.refs[idx - 1]
        if t == 1:                                      # SYMSXP
            name = self.item(); o = RObj("symbol", name); self.refs.append(o); return o
        if t == 9:                                      # CHARSXP
            n = self.i()
            return None if n == -1 else self.raw(n).decode("utf-8", "replace")
        if t in (2, 6, 17, 5):                          # LISTSXP, LANGSXP, DOTSXP, PROMSXP-like pairlists
            attrs = self.attrlist() if has_attr else {}
            tag = self.item() if has_tag else None
            car = self.item(); cdr = self.item()
            items = [(tag.value if isinstance(tag, RObj) and tag.type == "symbol" else None, car)]
            if isinstance(cdr, RObj) and cdr.type in ("pairlist", "lang"):
                items += cdr.value
            return RObj("lang" if t == 6 else "pairlist", items, attrs)
        if t == 4:                                      # ENVSXP
            locked = self.i(); o = RObj("env"); self.refs.append(o)
            o.value = dict(enclos=self.item(), frame=self.item(), hashtab=self.item()); a = self.item()
            return o
        if t in (10, 13):                               # LGLSXP, INTSXP
            n = self.i(); v = list(struct.unpack(f">{n}i", self.raw(4 * n)))
            v = [None if x == NA_INT else x for x in v]
            return RObj("lgl" if t == 10 else "int", v, self.attrlist() if has_attr else {})
        if t == 14:                                     # REALSXP
            n = self.i(); v = list(struct.unpack(f">{n}d", self.raw(8 * n)))
            return RObj("real", v, self.attrlist() if has_attr else {})
        if t == 16:                                     # STRSXP
            n = self.i(); v = [self.item() for _ in range(n)]
            return RObj("str", v, self.attrlist() if has_attr else {})
        if t in (19, 20):                               # VECSXP, EXPRSXP
            n = self.i(); v = [self.item() for _ in range(n)]
            return RObj("list", v, self.attrlist() if has_attr else {})
        if t == 24:                                     # RAWSXP
            n = self.i(); return RObj("raw", self.raw(n), self.attrlist() if has_attr else {})
        if t == 15:                                     # CPLXSXP
            n = self.i(); self.raw(16 * n); return RObj("cplx", None, self.attrlist() if has_attr else {})
        if t in (3,):                                   # CLOSXP
            attrs = self.attrlist() if has_attr else {}
            env = self.item(); formals = self.item(); body = self.item()
            return RObj("closure", None, attrs)
        if t in (8, 7):                                 # BUILTINSXP / SPECIALSXP
            n = self.i(); self.raw(n); return RObj("builtin")
        if t == 25:                                     # S4SXP
            return RObj("s4", None, self.attrlist() if has_attr else {})
        if t == 238:                                    # ALTREP_SXP
            info = self.item(); state = self.item(); attr = self.item()
            cls = info.value[0][1].value if info else None
            # compact sequences: state = c(n, start, step)
            if cls in ("compact_intseq", "compact_realseq"):
                n, start, step = [int(x) for x in state.value]
                v = [start + k * step for k in range(n)]
                return RObj("int" if cls == "compact_intseq" else "real", v, self.pl_to_dict(attr))
            if cls in ("wrap_real", "wrap_integer", "wrap_string", "wrap_logical"):
                inner = state.value[0]; inner.attrs.update(self.pl_to_dict(attr)); return inner
            if cls == "deferred_string":
                inner = state.value[0] if isinstance(state, RObj) and state.type == "pairlist" else state
                return RObj("str", [str(x) for x in inner.value], self.pl_to_dict(attr))
            raise ValueError(f"unsupported ALTREP class {cls}")
        if t == 247:                                    # BCREPDEF etc. are not expected
            raise ValueError("bytecode not supported")
        raise ValueError(f"unsupported SEXP type {t} at {self.p}")
    def pl_to_dict(self, pl):
        if pl is None: return {}
        return {k: v for k, v in pl.value}
    def attrlist(self):
        return self.pl_to_dict(self.item())

def read_rds(path):
    raw = open(path, "rb").read()
    if raw[:2] == b"\x1f\x8b": raw = gzip.decompress(raw)
    elif raw[:3] == b"BZh": raw = bz2.decompress(raw)
    elif raw[:6] == b"\xfd7zXZ\x00": raw = lzma.decompress(raw)
    r = Reader(raw); r.header(); return r.item()

def as_dataframe(o):
    """R data.frame -> (rownames, {colname: values})."""
    names = o.attr("names").value
    rn = o.attr("row.names")
    rows = rn.value if rn is not None else None
    if rn is not None and rn.type == "int" and len(rows) == 2 and rows[0] is None:
        rows = list(range(1, abs(rows[1]) + 1))
    return rows, {n: c.value for n, c in zip(names, o.value)}
