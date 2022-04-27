from subprocess import PIPE, Popen
import os, shutil

def run_cmd(cmd, verbose=False, output=False,error=False):
    stream=Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = stream.communicate()
    
    stdout=stdout.decode('utf-8')
    stderr=stderr.decode('utf-8')
    
    if stderr:
        print(stderr, flush=True)
    
    if verbose:
        print(stdout, flush=True)
        
        
    if output:
        return stdout
    if error:
        return stderr
    
def split_list(l,n=1000):
    i=0    
    chunk = l[i*n:(i+1)*n]
    while chunk:
        yield chunk
        i+=1
        chunk = l[i*n:(i+1)*n]
        
def get_attr(f,suffix):
    keys = []
    f.visit(lambda key : keys.append(f[key].attrs[suffix]) if suffix in f[key].attrs else None)
    
    return keys[0]        