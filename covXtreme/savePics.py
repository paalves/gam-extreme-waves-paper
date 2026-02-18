import matplotlib.pyplot as plt
import os

def savePics(Nme, width=11, height=6, typ='png'):
    """
    Save the current figure to a file with specified dimensions and format.

    Args:
        Nme (str): Filename (with or without extension).
        width (float, optional): Saved width of figure in inches. Defaults to 11.
        height (float, optional): Saved height of figure in inches. Defaults to 6 (9*2/3).
        typ (str, optional): Filetype ('png', 'jpg', 'pdf'). Defaults to 'png'.
    """
    
    # Get current figure
    fig = plt.gcf()
    
    # Set figure size in inches
    fig.set_size_inches(width, height)
    
    # Try to set font to Helvetica to match MATLAB code, fallback to sans-serif
    try:
        plt.rcParams['font.family'] = 'Helvetica'
    except:
        plt.rcParams['font.family'] = 'sans-serif'

    # Determine DPI and format specific settings
    dpi = 150
    if typ == 'pdf':
        dpi = 300
    elif typ == 'jpg' or typ == 'jpeg':
        typ = 'jpg' # normalization
        dpi = 150
    elif typ == 'png':
        dpi = 150
        
    # Ensure filename has the correct extension
    if not Nme.lower().endswith(f".{typ}"):
        filename = f"{Nme}.{typ}"
    else:
        filename = Nme

    # Check for directory existence (robustness)
    directory = os.path.dirname(filename)
    if directory and not os.path.exists(directory):
        os.makedirs(directory)

    # Save the figure
    # bbox_inches='tight' is used to mimic MATLAB's behavior of fitting 
    # the content within the paper size without excessive whitespace.
    plt.savefig(filename, format=typ, dpi=dpi, bbox_inches='tight')
    
    # Optional: Close plot to free memory if generating many plots
    # plt.close(fig)