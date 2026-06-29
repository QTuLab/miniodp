import sys
import subprocess
import time
import os
import logging
import signal

logging.basicConfig(format='[%(asctime)s %(levelname)s] %(message)s', stream=sys.stdout)
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def run_shell(cmd):
    p = subprocess.Popen(['/bin/bash', '-o', 'pipefail'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, preexec_fn=os.setsid)
    pid = p.pid
    pgid = os.getpgid(pid)
    log.info('run_shell_cmd: PID={}, PGID={}, CMD={}'.format(pid, pgid, cmd))
    t0 = time.perf_counter()
    stdout, stderr = p.communicate(cmd)
    rc = p.returncode
    t1 = time.perf_counter()
    err_str = 'PID={pid}, PGID={pgid}, RC={rc}, DURATION_SEC={dur:.1f}, STDERR={stde}, STDOUT={stdo}'.format(pid=pid, pgid=pgid, rc=rc, dur=t1 - t0, stde=stderr.strip(), stdo=stdout.strip())
    if rc:
        try:
            os.killpg(pgid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        raise Exception(err_str)
    log.info(err_str)
    return stdout.strip('\n')

def layout_check(dataDir, sra):
    srafile = f'{dataDir}/srafile/{sra}/{sra}.sra'
    lay = run_shell(f'fastq-dump -X 1 --split-spot -Z {srafile} | wc -l')
    layout = 'paired' if lay == '8' else 'single' if lay == '4' else 'unknown'
    return layout

def require_file(path, description):
    if not os.path.isfile(path):
        raise SystemExit(f"Missing {description}: {path}")

def require_dir(path, description):
    if not os.path.isdir(path):
        raise SystemExit(f"Missing {description}: {path}")
