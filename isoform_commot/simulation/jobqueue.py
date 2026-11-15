"""Run a job list with a concurrency limit read from a file before every start (plan pde_r0 section 11).

  python -m simulation.jobqueue --jobs JOBS.tsv --limit-file LIMIT --log-dir DIR

JOBS.tsv (simulation.benchmark jobs): name, done (file whose existence marks the job finished), command.
Every job runs under `nice -n 19` with 2 threads (OMP/OPENBLAS/MKL). Jobs whose done file exists are skipped.
The limit file holds one integer; editing it changes the concurrency without a restart (0 pauses new starts).
status.tsv in the log directory records start and end time, pid and return code of every job.
A restarted queue adopts jobs that an earlier queue started and that are still running (status.tsv, /proc):
they count towards the limit and are not started again; failed jobs (done file missing) are run again.
"""
import argparse
from datetime import datetime, timezone
import os
from pathlib import Path
import shlex
import subprocess
import time

import pandas as pd

THREADS = dict(OMP_NUM_THREADS='2', OPENBLAS_NUM_THREADS='2', MKL_NUM_THREADS='2')


def stamp():
    return datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M:%S')


def read_limit(path, default):
    try:
        return max(int(Path(path).read_text().strip()), 0)
    except (OSError, ValueError):
        return default


def alive(pid):
    try:
        with open(f'/proc/{pid}/cmdline', 'rb') as f:
            return b'simulation.benchmark' in f.read()
    except OSError:
        return False


def still_running(status_path):
    """{job name: pid} for jobs whose last status line is a start and whose process is still alive."""
    last = {}
    if Path(status_path).exists():
        for line in Path(status_path).read_text().splitlines():
            parts = line.split('\t')
            if len(parts) >= 4 and parts[1] != 'queue':
                last[parts[1]] = (parts[2], parts[3])
    return {name: int(pid) for name, (event, pid) in last.items() if event == 'start' and alive(int(pid))}


def main(args):
    jobs = pd.read_csv(args.jobs, sep='\t', dtype=str)
    log_dir = Path(args.log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    status = open(log_dir/'status.tsv', 'a', buffering=1)
    env = dict(os.environ, **THREADS)
    adopted = still_running(log_dir/'status.tsv')          # started by an earlier queue and still alive
    pending = [j for j in jobs.itertuples(index=False) if not Path(j.done).exists() and j.name not in adopted]
    running = {}
    status.write(f'{stamp()}\tqueue\tstart\t{os.getpid()}\t{len(pending)} pending of {len(jobs)}, '
                 f'{len(adopted)} adopted\n')
    limit = 0
    while pending or running or adopted:
        for name, (proc, start) in list(running.items()):
            code = proc.poll()
            if code is not None:
                status.write(f'{stamp()}\t{name}\tend\t{proc.pid}\t{code}\t{round(time.time() - start)}s\n')
                del running[name]
        for name, pid in list(adopted.items()):
            if not alive(pid):
                status.write(f'{stamp()}\t{name}\tend\t{pid}\t?\tadopted\n')
                del adopted[name]
        limit = read_limit(args.limit_file, limit)
        while pending and len(running) + len(adopted) < limit:
            job = pending.pop(0)
            if Path(job.done).exists():
                continue
            log = open(log_dir/(job.name.replace('/', '__') + '.log'), 'a')
            proc = subprocess.Popen(['nice', '-n', '19'] + shlex.split(job.command), stdout=log, stderr=subprocess.STDOUT,
                                    stdin=subprocess.DEVNULL, env=env, cwd=args.cwd)
            running[job.name] = (proc, time.time())
            status.write(f'{stamp()}\t{job.name}\tstart\t{proc.pid}\n')
        time.sleep(args.poll)
    status.write(f'{stamp()}\tqueue\tdone\t{os.getpid()}\n')


def parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--jobs', required=True)
    p.add_argument('--limit-file', required=True)
    p.add_argument('--log-dir', required=True)
    p.add_argument('--cwd', default='.')
    p.add_argument('--poll', type=float, default=5.0)
    return p


if __name__ == '__main__':
    main(parser().parse_args())
