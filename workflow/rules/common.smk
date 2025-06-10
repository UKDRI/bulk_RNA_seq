from pathlib import Path
import resource
import time


def adjust_ulimit():
    """ Adjust maximum allowed of number of file descriptor """
    soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    new_limit = min(hard, 65536)  # Set to 65536 or the maximum allowed
    if soft < new_limit:
        resource.setrlimit(resource.RLIMIT_NOFILE, (new_limit, hard))


def get_star_tempdir(wildcards, output):
    """ Create a temporary directory name based on sample and current time """
    return Path(output.bam).parent / ".".join([wildcards.sample, str(time.time())])