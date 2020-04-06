from utils.utils import make_valid_csv_file, make_valid_dir_string, make_dir
import queue
import time

SLEEPTIME = 0.001

def write_iter_of_dicts(iter_of_dicts, file_name, no_header=False):
    '''
    Write a csv from an iterable of dictionaries with entries as column names
    Dictionaries should only be one level deep, no nested dictionary

    Inputs:
        iter_of_dicts: an iterable of dictionaries to write
        file_name: name of file to write to
    kwargs:
        no_header: bool don't use a header line. Default=False
    Outputs: 
        string file name written to
    '''
    if len(iter_of_dicts) == 0:
        return 

    d = '/'.join(file_name.split('/')[:-1])
    d = make_valid_dir_string(d)
    make_dir(d)
    file_name = make_valid_csv_file(file_name)
    d = None 
    keys = [key for key in iter_of_dicts[0]]
    template = ''
    for _ in range(len(keys)):
        template += '{},'
    template = template[:-1] + '\n'
    with open(file_name, 'w') as o:
        if not no_header:
            o.write(template.format(*keys))
        for d in iter_of_dicts:
            vals = [i for _, i in d.items()]
            o.write(template.format(*vals))     

    return(file_name) 


def threaded_writer_interface(q: queue.Queue, flist: list) -> None:
    '''
    Threaded interface for the csv writer

    Inputs:
        q:      queue.Queue object for threadsafe operations
                For normal Queue operations, each entry should be of the form (writeable, filename)
                To end this thread, pass in (None, None)
        flist:  list to append written files to. Python lists are threadsafe by default
    Outputs:
        None
    '''
    while True:
        if not q.empty():
            writeable, filename = q.get()
            if writeable is None and filename is None: 
                break
            filename = write_iter_of_dicts(writeable, filename)
            flist.append(filename)
            q.task_done()
        else: 
            time.sleep(SLEEPTIME)


