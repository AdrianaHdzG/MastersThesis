
import matplotlib.pyplot as plt
import scipy
import sys


def plot_disp(x_coord, U, U_dir, beamX, model, model_el, modelNA, model_elNA, fig_num):
    fig, ax = plt.subplots()
    ax.plot(x_coord, U * 1000, 'rs', label='GF')  # 'rs' for red squares
    ax.plot(beamX, model * 1000, color='skyblue', label='ASREpy_ep')  # 'skyblue' for light blue
    ax.plot(beamX, model_el * 1000, 'b-', label='ASREpy_el')  # 'b-' for blue line
    ax.plot(beamX, modelNA * 1000, color='skyblue', marker='o', label='ASREpy_ep_NA')  # 'skyblue' with circles
    ax.plot(beamX, model_elNA * 1000, 'bo-', label='ASREpy_el_NA')  # 'bo-' for blue line with circles
    ax.set_xlabel('x coordinate')
    ax.set_ylabel(f'{U_dir} [mm]')
    ax.set_title(f'Fig. {fig_num}: Plot of X-coordinate [m] vs {U_dir} [mm]')
    ax.legend()  # Add legend
    plt.show()
    import scipy.io

    # Creating a dictionary with data
    data = {
        'x_coordinate': x_coord,
        'U_disp': U,
        'U_direction': U_dir,
        'beamX': beamX,
        'model': model,
        'model_el': model_el,
        'modelNA': modelNA,
        'model_elNA': model_elNA,
        'fig_num': fig_num
    }
    # Save the data to a .mat file
    scipy.io.savemat(f"Plot{fig_num}py.mat", data)


# Function to plot the forces and moments

def plot_strain(model, beamX):
    fig, ax = plt.subplots()

    ax.plot(beamX, model.axialForce[:len(beamX)], color='skyblue', label='Axial Force')
    ax.plot(beamX, model.shearForce[:len(beamX)], color='blue', label='Shear Force')

    # Plot moments for nodes 1-2 as they share element 1, moment 1 in node one and 2 in node two
    for i in range(len(beamX)):
        ax.plot([beamX[i], beamX[i]], [model.moment[0,i], model.moment[1,i]], color='k', label=f'Moment {i+1}' if i == 0 else "")

    ax.set_xlabel('x coordinate')
    ax.set_ylabel('N [N], S [N], M [Nm]')
    ax.set_title('Plot of X-coordinate [m] vs internal force')
    ax.legend() # Add legend

    plt.show()




'''
def plot_disp_veri(x_coord, U, U_dir, beamX, model, model_el, modelNA, model_elNA, fig_num, ValiModel):
    sys.path.append('FunctionScripts')
    data = scipy.io.loadmat(ValiModel)
    # Extract the array
    vali_model = data[ValiModel]
    print(vali_model)

    fig, ax = plt.subplots()
    ax.plot(x_coord, U * 1000, 'rs', label='GF')  # 'rs' for red squares
    ax.plot(beamX, model * 1000, color='skyblue', label='ASREpy_ep')  # 'skyblue' for light blue
    ax.plot(beamX, model_el * 1000, 'b-', label='ASREpy_el')  # 'b-' for blue line
    ax.plot(beamX, modelNA * 1000, color='skyblue', marker='o', label='ASREpy_ep_NA')  # 'skyblue' with circles
    ax.plot(beamX, model_elNA * 1000, 'bo-', label='ASREpy_el_NA')  # 'bo-' for blue line with circles

    ax.set_xlabel('x coordinate')
    ax.set_ylabel(f'{U_dir} [mm]')
    ax.set_title(f'Fig. {fig_num}: Plot of X-coordinate [m] vs {U_dir} [mm]')
    ax.legend()  # Add legend
    plt.show()
'''
