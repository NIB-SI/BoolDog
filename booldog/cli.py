import os, sys
import argparse


















def arguments():
    parser = argparse.ArgumentParser(
        description='CLI interface to BoolDoG. ',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    


    return parser








if __name__ == '__main__':
    parser = arguments()
    args = parser.parse_args()
