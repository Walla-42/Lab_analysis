# imports for necessary functions
import numpy as np
from scipy.integrate import simps

# import .csv file and find the index values of our area of interest
df = np.loadtxt('CorrChrom.csv', delimiter=',')
print(np.where(df == 97.8))
print(np.where(df == 128.6))

# ---------------------------------------------Trapezoidal Rule--------------------------------------------------------
# my formula
start = df[489, 0] # define the start and end values and set them to a variable
end = df[643, 0]
n = 643 - 489 # define the number of trapezoids. Maximum number is the difference in our indexes
h = (end - start) / n # this is our change in x value in the trapezoidal rule formula
sf = (df[489, 1] + df[643, 1]) # add together f(x) and f(x+1) which aren't multiplied by 2 in this formula
for i in range(1, n):
    index = int(489 + i) # For loop increases our index by 1 every loop then recalls the value and adds it to the
    # counter function
    sf += 2 * (df[index, 1])
sf *= h / 2 # Multiplies the delta x value to the counter function and divides by 2

print('Trapezoidal Integration: ' + str(sf))

# Check that our calculated value matches the Numpy Functions output
if round(sf, 10) == np.trapz(df[489:643 + 1, 1], df[489:643 + 1, 0]):
    print('True')
else:
    print('False')

# -----------------------------------------Simpsons method of integration-----------------------------------------------
# My function
a_index = 489  # Set index values a and b
b_index = 643
a = df[a_index - 1, 0]  # assign values to min and max of simpson formula f(x) and f(x+1)
b = df[b_index - 1, 0]
n1 = b_index - a_index  # assign n value or number of integrations to a variable
h1 = (b - a) / n1  # assign change of x to a variable
k = 0.0  # This is our counter variable. We need to set it to zero. This is especially important if this were a def
t = a + h1  # This is our second counter variable that will help us differentiate between the 2*f(x) and 4*f(x) values

for i in range(1, n1):
    if i % 2 == 0:  # Even Numbers will get the 2*f(x)
        k += 2 * df[a_index + i - 1, 1]
    else:  # Odd Numbers will get the 4*f(x)
        k += 4 * df[a_index + i - 1, 1]
    t += h1  # increases by h1 for every value in our formula from 1 to n1-1

result = (h1 / 3) * (df[a_index - 1, 1] + df[b_index - 1, 1] + k)  # Final Formula

print('Simpsons method integration: ' + str(result))

# -----------------------------------------------Scipy Simpsons Function-----------------------------------------------

print('Simpsons integration function: ' + str(simps(df[489:643, 1], df[489:643, 0])))  # Simpsons method is calculated
# using the scipy.simps formula to check our above function which was relatively accurate.

# ------------------------------------------------Numpy Trapz function-------------------------------------------------
trpz = np.trapz(df[489:643 + 1, 1], df[489:643 + 1, 0])
print('np.trap Funcion Integration: ' + str(trpz)) # We generated the np.trap formula to verify our above values
