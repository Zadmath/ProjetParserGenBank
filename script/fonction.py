#!/usr/bin/env python3
from Bio import Entrez
from lxml import etree
import pandas as pd
import os
import shutil
from Bio import SeqIO
import pickle
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
import random
import string
import signal
from ftplib import FTP
import re
import time
import os.path
import datetime
from threading import Thread
import functools