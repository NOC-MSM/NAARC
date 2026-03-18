conda activate pybdy

JVM="jre-17"
export JAVA_HOME=/usr/lib/jvm/$JVM/
cd /home/users/benbar/work/pybdy_examples/NAARC
 
run pybdy -s ./namelist_naarc.bdy
