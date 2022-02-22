#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 21 00:26:13 2021

@author: Ellie Pearce
"""


## Process ##

# 1. get mcfost output data
# 2. make source file
# 3. change data into astrochem input file
# 4. input file into astrochem and get abundances output
# 5. re input abundances into mcfost


### GET MCFOST OUTPUT ###

def get_mcfost_output(datadir):
    """

    Function take output of mcfost, all in one file, and converts it into
    required parameters for astrochem input file.

    """
    from astropy.io import fits
    from astropy import units
    import os
    from astropy import constants
    import numpy as np

    datadir = os.path.normpath(os.path.expanduser(datadir))

    ## Gas density ##
    hdu = fits.open(datadir + "/gas_density.fits.gz")
    gas_mass_density = (hdu[0].data * units.g *
                        units.cm**-3)   # g cm^-3
    gas_number_density = gas_mass_density / (2.0 * constants.u.to('g'))  # cm-3
    hdu.close()

    ## Dust temperature ##
    hdu = fits.open(datadir + "/Temperature.fits.gz")
    temperature = hdu[0].data * units.K
    hdu.close()

    ## UV field ##
    hdu = fits.open(datadir + "/UV_field.fits.gz")
    chi = hdu[0].data  # Habing units; convert into Draine units
    hdu.close()

    ## Grain mass ##
    hdu = fits.open(datadir + "/dust_mass_density.fits.gz")
    dust_mass_density = (hdu[0].data * units.g *
                         units.cm**-3)    # g cm^-3
    dust_gas_mass_ratio = dust_mass_density / gas_mass_density
    hdu.close()

    ## Grain size ##
    hdu = fits.open(datadir + "/grain_sizes.fits.gz")
    grain_sizes = hdu[0].data * units.um  # um
    hdu.close()
    hdu = fits.open(datadir + "/dust_particle_density.fits.gz")
    dust_number_density = (hdu[0].data *
                           units.m**-3)   # m^-3 per grain size bin
    hdu.close()
    average_grain_size = np.sqrt(sum(grain_sizes**2 *
                                     dust_number_density) /
                                 sum(dust_number_density)).to('um')

    # X-ray ionization (Bai & Goodman 2009) ##
    zeta1 = 6e-12 * units.s**-1  # s-1 (Tx = 3 keV)
    N1 = 1.5e21 * units.cm**-2   # cm
    zeta2 = 1e-15 * units.s**-1  # s-1
    N2 = 7e23 * units.cm**-2     # cm
    Lx29 = 5.                    # 10^29 erg s-1
    alpha = 0.4
    beta = 0.65
    hdu = fits.open(datadir + "/Column_density.fits.gz")
    column_mass_density_h = (hdu[0].data[0, :, :] *
                             units.g * units.cm**-2)   # g cm^-2, from the star
    column_density_h = column_mass_density_h / \
        (2.0 * constants.u.to('g'))  # cm-2
    column_mass_density_v = (hdu[0].data[1, :, :] *
                             units.g * units.cm**-2)   # g cm^-2, from the disk surface
    column_density_v = column_mass_density_v / \
        (2.0 * constants.u.to('g'))  # cm-2
    hdu.close()
    hdu = fits.open(datadir + "/grid.fits.gz")
    radius_au = hdu[0].data  # au
    hdu.close()
    zeta_xray = Lx29 / radius_au**2.2 * \
        (zeta1 * (np.exp(-(column_density_v/N1)**alpha)) +
         zeta2 * (np.exp(-(column_density_h/N2)**beta)))  # s-1 per H
    zeta_xray /= 2                                               # s-1 per H2

    # Cosmic-ray ionization (Bai & Goodman 2009) ##
    zeta0 = 1.3e-17 * units.s**-1  # "standard" value
    zeta_cr = zeta0 * \
        (np.exp(-(column_mass_density_v / (96 * units.g * units.cm**-2))))

    class Model:
        pass
    model = Model()
    model.gas_num_den = gas_number_density
    model.temperature = temperature
    model.chi = chi
    model.dust_gas_m_ratio = dust_gas_mass_ratio
    model.effective_grain_size = average_grain_size
    model.zeta_x = zeta_xray
    model.cosmic = zeta_cr
    model.density = gas_mass_density
    model.grain_abundance = dust_number_density

    return (model)


def source_file(model, index, source_filename):
    """


    Makes Astrochem source file.


    """

    f = open(source_filename, 'w')

    f.write('%s    0    %s     %s    %s' % (index, model.gas_num_den.ravel(
    ).value[index], model.temperature.ravel().value[index], model.temperature.ravel().value[index]))

    f.close()

    return source_filename


def make_astro_input(model, index, filename, abundances, output, sourcefile):
    """

    Take parameter model from get_mcfost_output and converts it into astrochem
    input init file.

    'model' is the model produced by get_mcfost_output

    'filename' is the name you want the outputted astrochem input file to have
    must have .ini at the end i.e. 'test.ini'

    'abundances' must be of format [[species name, initial abundance]]
    i.e CS --> abundance = [[CS, 1e-6]]

    'output' species of abundaces wanted in output i.e. ['CS,CO,S(+)']

    'sourcefile' needs to be in form 'source.mdl' from source_file()
    """
    # Getting values fomr model

    chi = model.chi.ravel()[index]
    cosmic = model.cosmic.ravel().value[index]
    grain_size = model.effective_grain_size.ravel().value[index]
    grain_gas_mass_ratio = model.dust_gas_m_ratio.ravel().value[index]
    grain_mass_density = model.density.value.ravel()[index]

    # Writes initial file headings
    name = filename+'.ini'

    f = open(name, 'w')

    f.write('[files]\n')
    f.write('source = %s\n' % (sourcefile))
    f.write('chem = osu2009.chm\n')
    f.write('[phys]\n')

    # Adding physical parameter values to input file

    f.write('%s = %s\n' % ('chi', chi))
    f.write('%s = %s\n' % ('cosmic', cosmic))
    f.write('%s = %s\n' % ('grain_size', grain_size))
    f.write('%s = %s\n' % ('grain_gas_mass_ratio', grain_gas_mass_ratio))
    f.write('%s = %s\n' % ('grain_mass_density', grain_mass_density))

    # Adding Solver parameters

    f.write('[solver]\n')
    f.write('abs_err = %s\n' % (10e-21))
    f.write('rel_err = %s\n' % (10e-6))

    # Initial abundances

    f.write('[abundances]\n')

    for i in range(len(abundances)):
        f.write('%s = %s\n' % (abundances[i][0], abundances[i][1]))

    # Output

    f.write('[output]\n')
    f.write('%s = %s\n' % ('abundances', output[0]))
    f.write('time_steps = 128\n')
    f.write('trace_routes = 1')

    # Close file

    f.close()

    return (name)


def run_astrochem(filename):
    """

    This function takes the input file made in make_astro_input and runs it
    in astrochem it then returns the the abundances file.

    'filename' is the name of the output from make_astro_input must have .ini 
    at the end i.e. 'test.ini'

    """
    import subprocess
    import os
    from os import environ

    env = dict(os.environ)
    env['DYLD_LIBRARY_PATH'] = '/Users/Ellie/opt/anaconda3/lib'

    # subprocess.run('DYLD_LIBRARY_PATH=/Users/Ellie/opt/anaconda3/lib; echo $DYLD_LIBRARY_PATH > crap.out', shell=True, capture_output=False, executable="/bin/zsh")
    # DYLD_LIBRARY_PATH=/Users/Ellie/opt/anaconda3/lib:/usr/local/lib
    output = subprocess.run("source ~/.zshrc;/usr/local/bin/astrochem " + filename, capture_output=False, shell= True, executable="/bin/zsh")
    if output.returncode:
        raise OSError(
            "astrochem did not run as expected, check astrochem's output")
    print("astrochem: Done")


def plot_abund(h5_file, abund_name, cell_num, abund_index=0):
    """


    Parameters
    ----------
    h5_file : NAME OF .H5 (ASTROCHEM OUTPUT) FILE
        MUST BE A 'STRING'

    abund_name : NAME OF ELEMENT YOU WNAT TO PLOT
        MUST BE A STRING

    cell_num : CELL NUMBER YOU ARE PLOTTING

    abund_index : INDEX OF NAME IF MORE THAN ONE ELEMENT
        BEING MODELLED


    Returns
    -------
    PLOT OF ABUNDANCE VS TIME WITH LOG-LOG SCALE

    """

    import h5py
    import matplotlib.pyplot as plt
    import numpy as np
    import subprocess

    # Import data
    f = h5py.File(h5_file, 'r')
    abund = np.array(f['Abundances'])
    time = np.array(f['TimeSteps'])
    species = np.array(f['Species'])
    name = '%s_%s' % (abund_name, cell_num)

    # Plot Data
    plt.plot(time, abund[abund_index][:])
    plt.title('%s Abundance: Cell %s (log-log scale)' % (abund_name, cell_num))
    plt.xlabel('Time')
    plt.ylabel('Abundance')
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig(name +'.pdf')
    plt.show()

    subprocess.run('mv %s.pdf abund_plots' % (name), shell=True)


def get_abund(data, ini_file, sourcefile, initial, chem, cell=[0], plot=0, plot_cell=0):
    
    """
    
    Parameters
    -----------
    data : mcfost output file in the form of a string
    
    ini_file : name you wnat ini file to have
    
    sourcefile : name you want source file to have
    
    initial : in the form: [['CS', 1e-6],['CO', 1e-4]]
    
    chem : chemical compound we want to find i.e. CS, inputted as 'CS'
    
    cell=[0] : is the cell we are finding the abundance for
    
    plot=0 : if you want the abundance to be plotted plot = 1, if you want it 
        to skip plotting plot = 0. Default is to not plot.

    plot_cell=0 : cell number you are plotting
    
    Returns
    --------
    ASTROCHEM OUTPUT AND OPTIONALLY A PLOT OF ABUNDANCES VS TIME (LOG-LOG SCALED)
        
    """

    import subprocess

    model = get_mcfost_output(data)

    for i in range(len(cell)):

        source = source_file(model, cell[i], sourcefile+'%s' % (cell[i]))

        file = make_astro_input(
            model, cell[i], ini_file + '%s' % (cell[i]), initial, chem, source)

        run_astrochem(file)

        if plot == 1:

            plot_abund(file, chem, cell, plot_cell)

        new_name = 'astrochem_output_%s.h5' % (cell[i])
        subprocess.run('mv astrochem_output.h5 ' + new_name, shell=True)

