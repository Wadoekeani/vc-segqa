"""Fetch a window of one z-plane at any pyramid level, without the whole slice.

section.py pulls a complete plane, which is fine at 45 um but not on the 2.4 um
volumes where a plane is hundreds of megabytes. Same trick - the volumes are
uncompressed, so a z-plane inside a chunk is a contiguous 16 KB block - applied
only to the chunks the window touches.
"""
import json, urllib.request, numpy as np
from concurrent.futures import ThreadPoolExecutor
from survey import BUCKET

CH = 128

def zarray(scroll, vol, level):
    return json.load(urllib.request.urlopen(
        f"{BUCKET}/{scroll}/volumes/{vol}/{level}/.zarray", timeout=30))

def plane_window(scroll, vol, level, z, cy, cx, half):
    """(image, y0, x0) covering [cy±half, cx±half] in that level's voxels."""
    shape = zarray(scroll, vol, level)["shape"]
    zc, zl = divmod(z, CH)
    y0c, y1c = max(0, int(cy - half)) // CH, min(shape[1] - 1, int(cy + half)) // CH
    x0c, x1c = max(0, int(cx - half)) // CH, min(shape[2] - 1, int(cx + half)) // CH
    jobs = [(y, x) for y in range(y0c, y1c + 1) for x in range(x0c, x1c + 1)]
    off = zl * CH * CH

    def get(t):
        y, x = t
        req = urllib.request.Request(f"{BUCKET}/{scroll}/volumes/{vol}/{level}/{zc}/{y}/{x}",
                                     headers={"Range": f"bytes={off}-{off + CH*CH - 1}"})
        try:
            return np.frombuffer(urllib.request.urlopen(req, timeout=90).read(),
                                 np.uint8).reshape(CH, CH)
        except Exception:
            return np.zeros((CH, CH), np.uint8)

    with ThreadPoolExecutor(16) as ex:
        parts = list(ex.map(get, jobs))
    img = np.zeros(((y1c - y0c + 1) * CH, (x1c - x0c + 1) * CH), np.uint8)
    for (y, x), p in zip(jobs, parts):
        img[(y - y0c) * CH:(y - y0c + 1) * CH, (x - x0c) * CH:(x - x0c + 1) * CH] = p
    return img, y0c * CH, x0c * CH, len(jobs) * CH * CH
