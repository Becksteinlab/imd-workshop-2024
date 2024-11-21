import matplotlib.pyplot as plt
from IPython.display import display, clear_output


def start_interactive_display():
    plt.ion()  # Turn on interactive mode


def stop_interactive_display():
    plt.ioff()


def update_display(obj):
    clear_output(wait=True)
    display(obj)
