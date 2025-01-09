"""
The module contains the function 'load_config' which loads the parameters from input file.
"""

def load_config(filename):
    """
    Loads the parameters from input file.

    Args:
        filename (str): Name of the input file.

    Returns:
        dict: The dictionary containing the loaded parameters.
    """
    config = {}
    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            key_value, _ = line.split('#', 1)

            key_value = key_value.strip()
            if '=' in key_value:
                key, value = key_value.split('=', 1)
                key = key.strip()
                value = value.strip()

                if value.isdigit():
                    config[key] = int(value)
                elif value.replace('.', '', 1).isdigit():
                    config[key] = float(value)
                else:
                    config[key] = value
    return config
