Vendor data is not redistributed. To score the measured inductor (scripts/vendor_chip.py):
download the 0201DS S-parameter archive from the Coilcraft 0201DS product page, extract 02DS-2N3.s2p, then run
    python scripts/coilcraft_to_csv.py path/to/02DS-2N3.s2p
which writes vendor/coilcraft_02DS-2N3_Zseries.csv here. reproduce.py does not need it.

MACOM diode S-parameters (README 'Diode model' numbers): download MA4AGP907_SPAR.zip and MA4AGFCP910_SPAR.zip from
https://cdn.macom.com/s-parameters/, unzip the .s2p files into vendor/macom/, then run  python scripts/macom_extract.py
