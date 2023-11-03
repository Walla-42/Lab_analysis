def WaterDensity(Temperature):
    # Define Variables used in Equation
    x = .99989
    x1 = 5.3322e-5
    x2 = -7.5899e-6
    x3 = 3.6719e-8
    # Equation used along with Temperature variable input with function
    density = (x + (x1 * Temperature) + (x2 * (Temperature ** 2)) + (x3 * (Temperature ** 3)))
    # Returns variable named density with an associated output value
    return density


def Bouyancy_Correction(Mass, T):
    dw = 8.00
    da = .0010
    x = .99989
    x1 = 5.3322e-5
    x2 = -7.5899e-6
    x3 = 3.6719e-8
    dH2O = (x + (x1 * T) + (x2 * (T ** 2)) + (x3 * (T ** 3)))
    Mass_Corrected = Mass * (1 - (da / dw)) / (1 - (da / dH2O))
    return Mass_Corrected


def EquivPoint(file,Burret_cal): #Returns the equivelence point for a potentiometric titration curve
    import numpy as np
    import matplotlib.pyplot as plt
    Vol = file[:, 0]* Burret_cal
    E = file[:, 1]
    # Define the derivatives of the parsed data
    dVol = np.diff(Vol)
    dE = np.diff(E)
    dxdy = dE / dVol
    # Define the average volume values
    AvgVol = (Vol[0:len(Vol) - 1] + Vol[1:len(Vol)]) / 2
    # Define index and value of the equivalence point from the derivative max and store as a float for manipulation
    peak = max(dxdy)
    devindex = np.where(dxdy == peak)[0]
    Equivpoint = AvgVol[devindex]
    Equivi = float(Equivpoint[0])
    # #You can even create a plot showing the data and the titration eq point with this code below:
    # fig, ax1 = plt.subplots()
    # ax2 = ax1.twinx()
    # ax1.plot(Vol, E, color='b', label='Titration')
    # ax2.plot(AvgVol, dxdy, color='r', label='Derivative')
    # ax1.set_xlabel('Electric Potential (V)')
    # ax1.set_ylabel('Volume Delivered (mL)')
    # ax2.set_ylabel('1st Derivative (mL/V)')
    # ax1.set_title('Equivalence Point Calculation')
    # fig.legend(loc=1)
    return Equivi

#This function is useful for when I need to input specific points. It will find the index of the point closes to my input value. 
def closest_index(data_frame, target):
    index = min(range(len(data_frame)), key=lambda i: abs(data_frame[i] - target))
    return index
