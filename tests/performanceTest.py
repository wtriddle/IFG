"""Script to test main program for its execution time"""

from main import main 
from datetime import datetime

now = datetime.now()
data = main()
print("exec_time: ", datetime.now() - now)