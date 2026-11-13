"""Download immutable official COMMOT source files and verify their checksums.

vendor/_unot.py           L1-UNOT core used by research.py (unchanged since v1).
vendor/commot_ot/         _unot.py + _cot.py as a tiny package, so the original
                          COMMOT collective OT (cot_combine_sparse) can be run as
                          a baseline without installing the full COMMOT package.
Upstream files are downloaded, not redistributed with this repository.
"""
import hashlib
from pathlib import Path
from urllib.request import urlopen

COMMIT = 'd117445bc07eaa19109c7609e97d9e35b26e99ca'
BASE = f'https://raw.githubusercontent.com/zcang/COMMOT/{COMMIT}/commot/_optimal_transport/'
URL = BASE + '_unot.py'
SHA256 = '4b50aca8e968277f04041532a78d4266b31ed4843937ea2e8e515aaae6527920'
COT_URL = BASE + '_cot.py'
COT_SHA256 = 'a19460aa08f9e18b9a9265ea8646f1e73284432f630bb604cbfbfc7acc987123'
VENDOR = Path(__file__).parent / 'vendor'
TARGET = VENDOR / '_unot.py'
PACKAGE = VENDOR / 'commot_ot'
FILES = {
    TARGET: (URL, SHA256),
    PACKAGE / '_unot.py': (URL, SHA256),
    PACKAGE / '_cot.py': (COT_URL, COT_SHA256),
}


def fetch():
    cache = {}
    for target, (url, sha) in FILES.items():
        if url not in cache:
            data = urlopen(url, timeout=60).read()
            if hashlib.sha256(data).hexdigest() != sha:
                raise RuntimeError(f'Upstream checksum mismatch: {url}')
            cache[url] = data
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(cache[url])
        print(target)
    # Package marker written locally; it only makes the relative import in _cot.py work.
    (PACKAGE / '__init__.py').write_text('"""Pinned COMMOT solver files; see fetch_solver.py."""\n')


if __name__ == '__main__':
    fetch()
