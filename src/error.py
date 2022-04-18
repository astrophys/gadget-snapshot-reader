import sys
def exit_with_error(string):
    """
    ARGS:
        string      : str to print then exit
    DESCRIPTION:
        Print string. Exit with value 1
    RETURN:
        N/A
    DEBUG:
        1. Tested, it worked
    FUTURE:
    """
    sys.stderr.write(string)
    sys.exit(1)


def warning(String):
    """
    ARGS:
        string      : string to print then exit
    DESCRIPTION:
        Print string. Exit with value 1
    RETURN:
        N/A
    DEBUG:
        1. Tested, it worked
    FUTURE:
    """
    sys.stderr.write(String)
