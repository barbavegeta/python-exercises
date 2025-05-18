import numpy as np
import sys

length = 0
names = ['Curie', 'Darwin', 'Turing']
for value in names:
    length = length + 1
    
print('There are', length, 'names in the list.')

# Examples

if '':
    print('empty string is true')
if 'word':
    print('word is true')
if []:
    print('empty list is true')
if [1, 2, 3]:
    print('non-empty list is true')
if 0:
    print('zero is true')
if 1:
    print('one is true')

############################################################################################

filenames = ['inflammation-01.csv',
         'myscript.py',
         'inflammation-02.csv',
         'small-01.csv',
         'small-02.csv']
large_files = []
small_files = []
other_files = []

for name in filenames:
    if name.startswith('inflammation'):
        large_files.append(name)
    elif name.startswith('small'):
        small_files.append(name)
    else:
        other_files.append(name)

print(other_files)
print(large_files)
print(small_files)

############################################################################################

name = 'name'
point = '*'

def fence(name, point):
    return point + name + point

fence(name, point)

############################################################################################

def add(a, b):
    print(a + b)

A = add(7, 3)
print(A)

############################################################################################

def add(a, b):
    return(a + b)

A = add(7, 3)
print(A)

############################################################################################

s = 'helium'

def outer(s):
    print(s[0] + s[-1])

B = outer(s)

def outer(s):
    return(s[0] + s[-1])

print(outer('helium'))

############################################################################################

def rescale(input_array):
    L = np.amin(input_array)
    H = np.amax(input_array)
    output_array = (input_array - L) / (H - L)
    return output_array

"""Takes an array as input, and returns a corresponding array scaled so
that 0 corresponds to the minimum and 1 to the maximum value of the input array.

Examples:
>>> rescale(numpy.arange(10.0))
array([ 0.        ,  0.11111111,  0.22222222,  0.33333333,  0.44444444,
       0.55555556,  0.66666667,  0.77777778,  0.88888889,  1.        ])
>>> rescale(numpy.linspace(0, 100, 5))
array([ 0.  ,  0.25,  0.5 ,  0.75,  1.  ])
"""

def rescale(input_array, low_val=0.0, high_val=1.0):
    """rescales input array values to lie between low_val and high_val"""
    L = np.amin(input_array)
    H = np.amax(input_array)
    intermed_array = (input_array - L) / (H - L)
    output_array = intermed_array * (high_val - low_val) + low_val
    return output_array

############################################################################################

def numbers(one, two=2, three, four=4):
    n = str(one) + str(two) + str(three) + str(four)
    return n

print(numbers(1, three=3))

############################################################################################

def favorite_ice_cream():
    ice_creams = [
        'chocolate',
        'vanilla',
        'strawberry'
    ]
    print(ice_creams[2])

favorite_ice_cream()

############################################################################################

def print_message(day):
    messages = [
        'Hello, world!',
        'Today is Tuesday!',
        'It is the middle of the week.',
        'Today is Donnerstag in German!',
        'Last day of the week!',
        'Hooray for the weekend!',
        'Aw, the weekend is almost over.'
    ]
    print(messages[day])

def print_sunday_message():
    print_message(6)

print_sunday_message()

############################################################################################

count = 0
for number in range(10):
    count = count + number
    print('The count is:', count)

############################################################################################

file_handle = open("my_file.txt", "r")
file_handle.read()

############################################################################################

def another_function():
    print('Syntax errors are annoying.')
    print('But at least Python tells us about them!')
    print('So they are usually not too hard to fix.')

still_function = another_function()

############################################################################################

message = ''

for number in range(10):
    # use a if the number is a multiple of 3, otherwise use b
    if (number % 3) == 0:
        message = message + 'a'
    else:
        message = message + 'b'

print(message)

############################################################################################

seasons = ['Spring', 'Summer', 'Fall', 'Winter']
print('My favorite season is', seasons[-1])

############################################################################################

numbers = [1.5, 2.3, 0.7, -0.001, 4.4]
total = 0.0
for num in numbers:
    assert num > 0.0, 'Data should only contain positive values'
    total += num
    print('total is:', total)

############################################################################################

def normalize_rectangle(rect):
    """Normalizes a rectangle so that it is at the origin and 1.0 units long on its longest axis.
    Input should be of the format (x0, y0, x1, y1).
    (x0, y0) and (x1, y1) define the lower left and upper right corners
    of the rectangle, respectively."""
    assert len(rect) == 4, 'Rectangles must contain 4 coordinates'
    x0, y0, x1, y1 = rect
    assert x0 < x1, 'Invalid X coordinates'
    assert y0 < y1, 'Invalid Y coordinates'
    dx = x1 - x0
    dy = y1 - y0
    if dx > dy:
        scaled = dx / dy
        upper_x, upper_y = 1.0, scaled
    else:
        scaled = dx / dy
        upper_x, upper_y = scaled, 1.0
    assert 0 < upper_x <= 1.0, 'Calculated upper X coordinate invalid'
    assert 0 < upper_y <= 1.0, 'Calculated upper Y coordinate invalid'
    return (0, 0, upper_x, upper_y)

print(normalize_rectangle( (0.0, 1.0, 2.0, 3.0) ))

############################################################################################

def range_overlap(ranges):
    """Return common overlap among a set of [left, right] ranges."""
    max_left = 0.0
    min_right = 1.0
    for (left, right) in ranges:
        max_left = max(max_left, left)
        min_right = min(min_right, right)
        if max_left >= min_right:
            return None
    return (max_left, min_right)

def test_range_overlap():
    assert range_overlap([ (0.0, 1.0), (5.0, 6.0) ]) == None
    assert range_overlap([ (0.0, 1.0), (1.0, 2.0) ]) == None
    assert range_overlap([ (0.0, 1.0) ]) == (0.0, 1.0)
    assert range_overlap([ (2.0, 3.0), (2.0, 4.0) ]) == (2.0, 3.0)
    assert range_overlap([ (0.0, 1.0), (0.0, 2.0), (-1.0, 1.0) ]) == (0.0, 1.0)
    assert range_overlap([]) == None

test_range_overlap()
#Traceback (most recent call last):
 # File "<stdin>", line 1, in <module>
  #File "<stdin>", line 5, in test_range_overlap
#AssertionError

#Traceback (most recent call last):
 # File "<stdin>", line 1, in <module>
  #File "<stdin>", line 2, in test_range_overlap
#AssertionError

############################################################################################

def get_total_cars(values):
    assert len(values) > 0
    for element in values:
        assert int(element)
    values = [int(element) for element in values]
    total = sum(values)
    assert total > 0
    return total

# The first assertion checks that the input sequence values is not empty. An empty sequence such as [] will make it fail.
# The second assertion checks that each value in the list can be turned into an integer. Input such as [1, 2, 'c', 3] will make it fail.
# The third assertion checks that the total of the list is greater than 0. Input such as [-10, 2, 3] will make it fail.

############################################################################################

patients = [[70, 1.8], [80, 1.9], [150, 1.7]]

def calculate_bmi(weight, height):
    print('weight:', weight, 'height:', height)
    return weight / (height ** 2)

for patient in patients:
    weight, height = patient
    bmi = calculate_bmi(weight, height)
    print("Patient's BMI (weight: %f, height: %f) is: %f" % (weight, height, bmi))
    
############################################################################################

print("sys version is", sys.version)
print("sys.argv is", sys.argv)

############################################################################################

################# FUNCTION

def main():
    script = sys.argv[0]
    action = sys.argv[1]
    filenames = sys.argv[2:]

    for filename in filenames:
        data = np.loadtxt(filename, delimiter=',')

        if action == '--min':
            values = np.amin(data, axis=1)
        elif action == '--mean':
            values = np.mean(data, axis=1)
        elif action == '--max':
            values = np.amax(data, axis=1)

        for val in values:
            print(val)

if __name__ == '__main__':
    main()

################# ASSERT


def main():
    script = sys.argv[0]
    action = sys.argv[1]
    if action not in ['--min', '--mean', '--max']: # if no action given
        action = '--mean'    # set a default action, that being mean
    filenames = sys.argv[2:]
    assert action in ['--min', '--mean', '--max'], \
           'Action is not one of --min, --mean, or --max: ' + action
    if len(filenames) == 0:
        process(sys.stdin, action)
    else:
        for filename in filenames:
            process(filename, action)

def process(filename, action):
    data = np.loadtxt(filename, delimiter=',')

    if action == '--min':
        values = np.amin(data, axis=1)
    elif action == '--mean':
        values = np.mean(data, axis=1)
    elif action == '--max':
        values = np.amax(data, axis=1)

    for val in values:
        print(val)

if __name__ == '__main__':
    main()

############################################################################################

def main():
    assert len(sys.argv) == 4, 'Need exactly 3 arguments'

    operator  = sys.argv[1]
    assert operator in ['--add', '--subtract', '--multiply', '--divide'], \
           'Action is not one of --add, --subtract, --multiply, or --divide: '
    try:
        operand1, operand2 = float(sys.argv[2]), float(sys.argv[3])
    except ValueError:
        print("Invalid input")
        return
    
    process(operand1, operator, operand2)

def process(operand1, operator, operand2):
    if operator == '--add':
        values = operand1 + operand2
    elif operator == '---subtract':
        values = operand1 - operand2
    elif operator == '--multiply':
        values = operand1 * operand2
    elif operator == '--divide':
        values = operand1 / operand2
    print(values)

main()

#############4

# this is code/readings_08.py

def main():
    script = sys.argv[0]
    if len(sys.argv) == 1:  # no arguments, so print help message
        print("Usage: python readings_08.py action filenames\n"
              "Action:\n"
              "    Must be one of --min, --mean, or --max.\n"
              "Filenames:\n"
              "    If blank, input is taken from standard input (stdin).\n"
              "    Otherwise, each filename in the list of arguments is processed in turn.")
        return

    action = sys.argv[1]
    filenames = sys.argv[2:]
    assert action in ['--min', '--mean', '--max'], (
        'Action is not one of --min, --mean, or --max: ' + action)
    if len(filenames) == 0:
        process(sys.stdin, action)
    else:
        for filename in filenames:
            process(filename, action)

def process(filename, action):
    data = np.loadtxt(filename, delimiter=',')

    if action == '--min':
        values = np.amin(data, axis=1)
    elif action == '--mean':
        values = np.mean(data, axis=1)
    elif action == '--max':
        values = np.amax(data, axis=1)

    for val in values:
        print(val)

if __name__ == '__main__':
    main()

#########################################################################################

def main():
    script = sys.argv[0]
    filenames = sys.argv[1:]
    if len(filenames) <=1: #nothing to check
        print('Only 1 file specified on input')
    else:
        nrow0, ncol0 = row_col_count(filenames[0])
        print('First file %s: %d rows and %d columns' % (filenames[0], nrow0, ncol0))
        for filename in filenames[1:]:
            nrow, ncol = row_col_count(filename)
            if nrow != nrow0 or ncol != ncol0:
                print('File %s does not check: %d rows and %d columns' % (filename, nrow, ncol))
            else:
                print('File %s checks' % filename)
        return

def row_col_count(filename): # .shape is what actually counts the number of rows and columns
    try:
        nrow, ncol = np.loadtxt(filename, delimiter=',').shape
    except ValueError:
        # 'ValueError' error is raised when numpy encounters lines that
        # have different number of data elements in them than the rest of the lines,
        # or when lines have non-numeric elements
        nrow, ncol = (0, 0)
    return nrow, ncol

main()

#########################################################################################


def main():
    """print each input filename and the number of lines in it,
       and print the sum of the number of lines"""
    filenames = sys.argv[1:]
    sum_nlines = 0 #initialize counting variable

    if len(filenames) == 0: # no filenames, just stdin
        sum_nlines = count_file_like(sys.stdin)
        print('stdin: %d' % sum_nlines)
    else:
        for filename in filenames:
            nlines = count_file(filename)
            print('%s %d' % (filename, nlines))
            sum_nlines += nlines
        print('total: %d' % sum_nlines)

def count_file(filename):
    """count the number of lines in a file"""
    f = open(filename, 'r')
    nlines = len(f.readlines())
    f.close()
    return(nlines)

def count_file_like(file_like):
    """count the number of lines in a file-like object (eg stdin)"""
    n = 0
    for line in file_like:
        n = n+1
    return n

main()