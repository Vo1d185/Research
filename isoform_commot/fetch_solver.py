"""Download one immutable official COMMOT source file and verify its checksum."""
import hashlib
from pathlib import Path
from urllib.request import urlopen

URL = 'https://raw.githubusercontent.com/zcang/COMMOT/d117445bc07eaa19109c7609e97d9e35b26e99ca/commot/_optimal_transport/_unot.py'
SHA256 = '4b50aca8e968277f04041532a78d4266b31ed4843937ea2e8e515aaae6527920'
TARGET = Path(__file__).parent / 'vendor' / '_unot.py'

if __name__ == '__main__':
    data = urlopen(URL, timeout=60).read()
    if hashlib.sha256(data).hexdigest() != SHA256:
        raise RuntimeError('Upstream solver checksum mismatch')
    TARGET.parent.mkdir(parents=True, exist_ok=True)
    TARGET.write_bytes(data)
    print(TARGET)
