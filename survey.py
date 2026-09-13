"""Which published scrolls can be QA-scanned? Needs segments + a masked volume."""
import sys, urllib.request, urllib.parse, xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor

BUCKET = "https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com"
NS = "{http://s3.amazonaws.com/doc/2006-03-01/}"

def list_dir(prefix):
    url = f"{BUCKET}/?list-type=2&delimiter=/&max-keys=1000&prefix={urllib.parse.quote(prefix)}"
    root = ET.fromstring(urllib.request.urlopen(url, timeout=60).read())
    return ([p.findtext(NS+"Prefix") for p in root.findall(NS+"CommonPrefixes")],
            [c.findtext(NS+"Key") for c in root.findall(NS+"Contents")])

def probe(scroll):
    name = scroll.rstrip("/").split("/")[-1]
    try:
        segs, _ = list_dir(scroll + "segments/")
        vols, _ = list_dir(scroll + "volumes/")
    except Exception as e:
        return name, 0, [], f"error {e}"
    masked = [v.rstrip("/").split("/")[-1] for v in vols if "masked" in v]
    return name, len(segs), masked, ""

if __name__ == "__main__":
    scrolls, _ = list_dir("")
    print(f"{len(scrolls)} top-level objects\n")
    with ThreadPoolExecutor(12) as ex:
        res = list(ex.map(probe, scrolls))
    usable = []
    print(f"{'scroll':22s} {'segs':>5}  masked volumes")
    for name, ns, masked, err in sorted(res, key=lambda r: -r[1]):
        if not ns and not masked: continue
        print(f"{name:22s} {ns:>5}  {', '.join(m[:46] for m in masked[:2]) or '-'}")
        if ns and masked: usable.append(name)
    print(f"\nscrolls with BOTH segments and a masked volume: {len(usable)}")
    print(" ".join(usable))
