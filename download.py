import requests
import sys
pyversion = sys.version_info.major
import urllib.request
import shutil
import tarfile as tar
import os
import glob
def download(url,filename):
    print('Downloading ' + url )
    if pyversion == 2:
        r = urllib2.urlopen(url).read()
        f = open(filename, 'w')   # write it to the right filename
        f.write(r)
        f.close()
    else:
        urllib.request.urlretrieve(url, filename)
    print("File download complete")
    return filename

#define EOS and pair of Masses to be downloaded (this includes every EOS and pair of masses)
EOS=['15H','125H','H','HB','B','SFHo']
MASS=['135_135','125_146','125_125','121_151','118_155','117_156','116_158','112_140','107_146']

#create folder data to save the downloaded files
if os.path.exists('data'):
    pass
else:
    os.mkdir('data')

#download and save in the appropriate folder with the appropriate name
for eos in EOS:
    for mas in MASS:
        for w in ['00155','0015']:
            for d in ['800','540']:
                url='http://www2.yukawa.kyoto-u.ac.jp/~nr_kyoto/SACRA_PUB/'+eos+'_'+mas+'_'+w+'_182_135/h_'+d+'.0_l2m2'
                request = requests.get(url)
                if request.status_code == 200:
                    download(url,'data/'+eos+'_'+mas)
                    break
