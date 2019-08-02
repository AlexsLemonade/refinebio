#!/usr/bin/env python3
# `buffer.py` reads stdin into an internal buffer and feeds that to stdout.
# It's kind of like a pipe, but the buffer has your system RAM as its capacity.

import queue
import sys
import threading
import time


queue = queue.Queue()
finished_reading = False
READ_SIZE_BYTES = 1024


def read_thread():
    global finished_reading
    while not finished_reading:
        output = sys.stdin.buffer.read(READ_SIZE_BYTES)
        if len(output) == 0:
            finished_reading = True
            sys.stdin.close()
        else:
            queue.put(output)


def write_thread():
    global finished_reading

    # We always want this loop to run once after finished_reading is true
    terminate = False
    while not terminate:
        terminate = finished_reading
        while not queue.empty():
            sys.stdout.buffer.write(queue.get())
        time.sleep(.05)

    sys.stdout.close()


if __name__ == "__main__":
    reading = threading.Thread(target=read_thread)
    reading.start()

    writing = threading.Thread(target=write_thread)
    writing.start()
