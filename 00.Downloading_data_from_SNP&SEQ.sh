#https://ngisweden.scilifelab.se/resources/data-delivery-dds/
module load bioinfo-tools
module load dds-cli
dds
#Here you add the username and password that you created when you received the registration e-mail. You will shortly receive an e-mail with a one-time authentication code (sent from services-noreply@scilifelab.se).
dds auth login
DDS username: fpertille
DDS password: same as mp inverted capital letters

#list projects
dds ls

dds ls -p snpseq00103 --tree
dds data get -p snpseq00103 -a --verify-checksum -d /proj/gbs_medip/GEroNIMO/raw

dds data get -p snpseq00168 -a --verify-checksum -d /proj/gbs_medip/uppstore2017266/private/raw

dds auth logout

snpseq00103
