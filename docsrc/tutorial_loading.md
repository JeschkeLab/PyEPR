# Loading and Analysing Data

The PyEPR package normally saves data in the HDF5 format, which is a widely used file format for storing large amounts of numerical data. Files which are saved in this format typically have the file extension `.h5` or `.hdf5`, and if created by PyEPR will have the metadata stored in an easily accessible format for PyEPR to read. Nonetheless, PyEPR can also load data from other file formats, such as text files (`.txt`, `.csv`), Bruker Elexsys files (`.dta`, `.dsc`).

Data is loaded into the xarray `Dataarray` object, which is a powerful data structure for handling multi-dimensional arrays with labeled axes and coordinates. This allows for easy manipulation, analysis, and visualization of the data.

## Loading Data

```python
import pyepr as epr
# Load data from an HDF5 file
data = epr.eprload('path/to/datafile.h5')
# Load data from a text file
data_txt = epr.eprload('path/to/datafile.txt')
```

## Viewing Data
Once the data is loaded, you can view its contents and metadata using standard xarray methods.

```python
# Print the data array
print(data)
# Access metadata
print(data.attrs)
```
The data can also be quickly visualized using the built-in plotting functions.

```python
# Plot the data
data.real.plot(label='Re')
data.imag.plot(label='Im')
```

Xarray has been extended to include some convient EPR methods, such as `correctphase`.

```python
# Correct the phase of the EPR signal
data_corrected = data.epr.correctphase()
data_corrected.plot(label='Corrected Re')
```

## Saving Data
After processing and analyzing the data, you can save it back to an HDF5 file for future use.
```python
# Save the processed data to an HDF5 file
data_corrected.epr.save('path/to/processed_datafile')
```

## PyEPR Automated Data Analysis
Since in EPR we often have to perform standard data analysis routines, PyEPR includes a set of automated data analysis functions that can be applied to the loaded data. These functions are designed to streamline the process of extracting the necessary information to move onto the next experiment. 
A list of currently implemented automated data analysis functions can be found in the [API documentation](API_docs.rst).

### Field Sweep Analysis

```python
fieldsweep_data = epr.eprload('path/to/fieldsweep_datafile.h5')
fieldsweep = epr.FieldSweepAnalysis(fieldsweep_data)
# Determine the gyromagnetic ratio
gyro = fieldsweep.calc_gyro()
# Plot the results
fieldsweep.plot()
```

### Relaxation Analysis

```{eval-rst}
.. tab-set::

    .. tab-item:: Repetition Time T1 Recovery

        .. code-block:: python

            t1_data = epr.eprload('path/to/t1_datafile.h5')
            t1_analysis = epr.ReptimeAnalysis(t1_data)
            # Calculate T1 relaxation time
            t1_result = t1_analysis.fit()
            # Plot the results
            t1_analysis.plot()

    .. tab-item:: Hahn Echo Tm Relaxation
        
        .. code-block:: python

            tm_data = epr.eprload('path/to/tm_datafile.h5')
            t2_analysis = epr.HahnEchoRelaxationAnalysis(t2_data)
            # Calculate T2 relaxation time
            t2_result = t2_analysis.fit()
            # Plot the results
            t2_analysis.plot()

```

### Resonator Profile Analysis
