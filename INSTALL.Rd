#############################################################
# This plugin requires Singularity deplopyed on TSS.
#
# We can test singularity if installed by typing "singularity"
#############################################################

$ singularity

#############################################################
# If not, we can install singularity easily:
#############################################################

# for ubuntu 14.04:

wget http://updates.itw/updates/software/trusty-dev/ion-singularity_3.5.3-2_amd64.deb
sudo dpkg -i ion-singularity_3.5.3-2_amd64.deb

# for ubuntu 18.04:
wget http://updates.itw/updates/software/bionic-dev/ion-singularity_3.5.3-2_amd64.deb
sudo dpkg -i ion-singularity_3.5.3-2_amd64.deb

# then test it:

$ singularity
