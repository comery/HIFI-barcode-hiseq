import sys
import numpy as np

def progress_bar(num_cur, total):
    ratio = float(num_cur) / total
    percentage = int(ratio * 100)
    r = '\r\n[%s%s]%d%%' % (">"*percentage, " "*(100-percentage), percentage )
    sys.stdout.write(r)
    sys.stdout.flush()

with open('hifi-hiseq.py', 'r') as f:
    lines = f.readlines()

print('Total line numbers are: {}'.format(len(lines)))
cur_ = 1
total_ = len(lines)
for line in lines:
    if cur_ % 10 == 0 or cur_ == total_:
        # your processing code here #
        progress_bar(cur_, total_)
    cur_ += 1
