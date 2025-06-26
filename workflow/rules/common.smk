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

def get_star_memory_gb(wildcards):
    """ Divide the value provided to limit star memory by a billion"""
    return config["star_ram_limit"] / 1000000000
    
def get_deseq_output_dir(wildcards, output):
    """ Create a temporary directory name based on sample and current time """
    return Path(output.deg_results[0]).parent
