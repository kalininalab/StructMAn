import ray
import os
from structman import settings
import subprocess
import time
import sys
import traceback

def ray_init(config, n_try = 0, redis_mem_default = False, overwrite_logging_level = None, total_memory_quantile = 0.5, debug = False):
    if ray.is_initialized():
        return
    
    if not config.python_env_expanded_for_ray:
        os.environ["PYTHONPATH"] = f'{settings.ROOT_DIR}:{os.environ.get("PYTHONPATH", "")}'
        os.environ["PYTHONPATH"] = f'{settings.LIB_DIR}:{os.environ.get("PYTHONPATH", "")}'
        os.environ["PYTHONPATH"] = f'{settings.RINERATOR_DIR}:{os.environ.get("PYTHONPATH", "")}'
        os.environ["PYTHONPATH"] = f'{settings.OUTPUT_DIR}:{os.environ.get("PYTHONPATH", "")}'
        
        config.python_env_expanded_for_ray = True

    logging_level = 20
    if config.verbosity <= 1:
        logging_level = 0

    if overwrite_logging_level is not None:
        logging_level = overwrite_logging_level

    if config.verbosity <= 3 and not debug:
        log_to_driver = False
    else:
        log_to_driver = True
        print(f'Setting log_to_driver in ray.init to True, local mode: {config.ray_local_mode}')

    redis_mem = ((total_memory_quantile)/10.0) * config.gigs_of_ram * 1024 * 1024 * 1024
 
    mem = ((total_memory_quantile*0.9)/10.0) * config.gigs_of_ram * 1024 * 1024 * 1024

    if config.verbosity >= 2:
        print(f'Call of ray.init: num_cpus = {config.proc_n}, object_store_memory = {redis_mem}')

    errs = '1'
    loops = 0
    while len(errs) > 0:
        try:
            p = subprocess.Popen(['ray','stop','--force'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            outs, errs = p.communicate()
        except:
            [e, f, g] = sys.exc_info()
            g = traceback.format_exc()
            config.errorlog.add_warning(f'Ray stop failed:\n{e}\n{f}\n{g}')
            outs = b''
            errs = b''

        loops += 1

        if config.verbosity >= 3:
            print(f'Returns of ray stop:\n{outs}\n{errs}')
        if outs.decode('ascii').count('Stopped') > 0:
            time.sleep(5)
        if loops > 1:
            time.sleep(10**loops)
        if loops >= 4:
            break
    try:
        if not redis_mem_default:
            ray.init(num_cpus=config.proc_n, include_dashboard=False, ignore_reinit_error=True, logging_level = logging_level,
                log_to_driver = log_to_driver, local_mode = config.ray_local_mode, object_store_memory = redis_mem, _memory = mem)
        else:
            ray.init(num_cpus=config.proc_n, include_dashboard=False, ignore_reinit_error=True, logging_level = logging_level,
                log_to_driver = log_to_driver, local_mode = config.ray_local_mode)
    except:
        [e, f, g] = sys.exc_info()
        g = traceback.format_exc()
        if n_try == 3:
            config.errorlog.add_error(f'Ray init failed:\n{e}\n{f}\n{g}')
            return

        config.errorlog.add_warning(f'Ray init failed, retry {n_try}... \n{e}\n{f}\n{g}')
        
        try:
            p = subprocess.Popen(['ray','stop','--force'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            outs, errs = p.communicate()
        except:
            [e, f, g] = sys.exc_info()
            g = traceback.format_exc()
            config.errorlog.add_warning(f'Ray stop failed:\n{e}\n{f}\n{g}')
            outs = None
            errs = ''

        if config.verbosity >= 3:
            print(f'Returns of ray stop:\n{outs}\n{errs}')
        time.sleep(10*n_try)
        if n_try <= 4:
           ray_init(config, n_try = (n_try + 1), total_memory_quantile = total_memory_quantile - 0.05)
        else:
           ray_init(config, n_try = (n_try + 1), redis_mem_default = True)

def ray_hack():
    # hack proposed by the devs of ray to prevent too many processes being spawned
    return
    #resources = ray.ray.get_resource_ids()
    #cpus = [v[0] for v in resources['CPU']]
    #psutil.Process().cpu_affinity(cpus)

def monitor_ray_store():
    p = subprocess.Popen(['ray','memory'])

    page, err = p.communicate()
    print(page)

    return
