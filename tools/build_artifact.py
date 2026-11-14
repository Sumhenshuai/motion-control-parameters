"""Build a single-page copy of app/index.html for the claude.ai artifact viewer.

The artifact viewer wraps the page in its own document, blocks service workers and
only allows scripts from a few CDNs, so this copy drops the document wrapper and
PWA tags and loads MathJax from jsDelivr instead of the bundled file.

Usage: python3 tools/build_artifact.py OUTPUT.html
"""
import sys
from pathlib import Path

src = Path(__file__).resolve().parent.parent / "app" / "index.html"
out = Path(sys.argv[1])
drop_exact = {
    '<!doctype html>', '<html lang="zh-CN">', '<head>', '<meta charset="utf-8">',
    '<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">',
    '</head>', '<body>', '</body>', '</html>',
}
drop_prefix = ('<link rel="manifest"', '<meta name="apple-', '<meta name="mobile-web-app', '<meta name="theme-color"',
               '<link rel="apple-touch-icon"', '<link rel="icon"')
lines = []
for line in src.read_text(encoding="utf-8").split("\n"):
    t = line.strip()
    if t in drop_exact or t.startswith(drop_prefix):
        continue
    if 'src="vendor/mathjax/tex-svg-full.js"' in t:
        line = '<script defer src="https://cdn.jsdelivr.net/npm/mathjax@3.2.2/es5/tex-svg-full.js"></script>'
    lines.append(line)
out.write_text("\n".join(lines), encoding="utf-8")
print(f"wrote {out}")
