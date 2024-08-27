#!/usr/bin/perl
#
# FOR N-subdomain motion, based on RMSD, Notice RMSD domain.
# Thus, use nprint = 200, to make 50ps for each state, total 100 ps.
#
# Modified for running on our node, not HPC
# ******  IMPORTANT  *****
# ******  IMPORTANT  *****
# ******  IMPORTANT  *****
#
# ****** A bug, $dummy = $nframe[1] - $chstp; is fixed as:
# ******        $dummy = $nframe[1] +1 - $chstp;
# ******	$nframe[2]+=($chstp)-1;
# ****** See details in Green Node
#
# ******  IMPORTANT  *****
# ******  IMPORTANT  *****
# ******  IMPORTANT  *****
#
# ****************** SECTION 1 ************************
# This section is specific to the current job
# The user is required to change the fields in this section
#
# Before run perl, need to prepare following:
# 1. path.
# 2. Build mass manually.
#    Use psf file, Excel. What is 2nd column for?
# 3. Modify tpsbarefile (input for NAMD and Amber-ptraj) accordingly.
# 4. For TMD raw data, prep _out and _bdd by running each fortran codes manually.
# 5. Put *.psf into ptraj.src
# 6. At # Applying velocity perturbtion #, put *.psf in grep section.
# 7. For this test.run 30ps, then require 30,000 steps with 1fs timestep.
#    Set $nprint for every 500 steps, then nframe[*]=60. Change nframe def at STEP 4.
#    This process may not be affected by initial nframe (at gen0) [need to confirm].
# 7A.Change new periodicity in tpsbarefile.
# 8. Change chstp = chstp X *
# 9. Prepare initial fake ene file for nframe computation at gen0.
#    (to mimic 0 output, total # of ENERGY: should be one more than nframe).
#10. Prepare *.dat6 file from initial ptraj (gen0) computation.
#11. Run tps_bol_ck manually to generate _out file.
#12. Run tps_bol_ck manually to generate _bdd file.
#
# Type of run TPS or BOLAS? set variable=1 for TPS, 2 for BOLAS
$tps_or_bolas=1;
#
# To terminate a run anytime type "cp end generation_status"
# at the command line and the program will exit after the current
# iteration
#
# jobname if using a queing system to submit jobs
$jobname="TPS_test";
#
# iteration file names (5 characters)
$job_file="gn_";
#
# Execution directory
$execdir="/scratch/sy3757/SimRNA_run/77I_mut36PSM/traj/01/TMD_77R/TPS_sampling/TPS_2500.cp-2";

# Path to directory which holds the tps_xxx_ck and
# bead_xxx_ck executables
$spiderexec="/scratch/sy3757/SimRNA_run/77I_mut36PSM/traj/01/TMD_77R/TPS_sampling/TPS_2500.cp-2";
# Define TOPOLIGY and extended files
$topology="step3_input.parm7";
# Define reference pdb file
#$reference ="step1_pdbreader.77nt.3_3.aligned.tmd.pdb";
#
# the executables are generated from fortran90 programs
# to generate executables, issue the following commands
# cd  $spiderexec
# f90 -o $dummy1 $dummy2
# f90 -o $dummy3 $dummy4
# where $dummy1=$spiderexec . "/tps_bol_ck"
#       $dummy2=$spiderexec . "/tps_bol_ck" . ".f";
#       $dummy3=$spiderexec . "/bead_bol_ck"
#       $dummy4=$spiderexec . "/bead_bol_ck" . ".f";

system "gfortran -o tps_bol_ck tps_bol_ck.f";
system "gfortran -o bead_bol_ck bead_bol_ck.f";
system "gfortran -o perturb perturb.f";
system "gfortran -o random random.f";
#system "gfortran -o convert convert.f";
#system "gfortran -o convert-vel convert-vel.f";

#
$shootexec=$spiderexec . "/tps_bol_ck";
$beadexec=$spiderexec . "/bead_bol_ck";
$goperturb=$spiderexec . "/perturb";
$random=$spiderexec . "/random";
#$goconvert=$spiderexec . "/convert";
#$goconvertvel=$spiderexec . "/convert-vel";
#

# command if using queing system (see section)
$communication="sbatch";
#
# Number of processors per job
$nproc=1;
# Number of processors are stored in file nproc. The value in
# nproc can be updated any time; the scripts reads the file nproc
# and reassigns the number of processors to be used

$countfile="nproc";
open (count4,">$countfile");
print count4 $nproc,"\n";
close(count4);

# 1 SHIFT for every FRQShift SHOOTS chosen at random. Integer value please.
$FRQshift=2;
#
# nsteps is half the trajectory length, i.e., (total number of timesteps)/2  # 10ns
$nsteps=1000000;
#
# time step in femto seconds
$timestep=0.001;
#
# template charmm input files to harvest forward and reverse trajectories
# during shooting/ shifting moves
$tpsbarefile="tps_amber";
#$tpsbarefile2="cpptraj_all";
$nopara=3;
#
# displacement step for shooting moves in units ...
#$dis_step=0.005; (Original)
#$dis_step=0.005;
$dis_step=0.1;
#
# Factor by which displacement step is altered if the acceptance rate differs
# from desired value. This only applies to TPS runs. For BOLAS, the
# acceptance rate is near unity, thus dis_fact is reset to 1.0 and
# $acc_desire is ignored
$dis_fact=1.05;
#
# desired value of acceptance rate (decimal)
$acc_desire=0.33;
#
# set the maximum value of shift as a fraction of total traj. length
# for example 0.15 implies 15% of trajectory length
$maximum_delta_shift=0.15;
#
# initialize seed for random generator (negative)
$irand=-18512;
#
# charmm prints and saves configurations every nprint steps in the
# trajectory files and the output file

#    Temporary # 1000 frames each traj
$nprint=2000;
#$nprint=10;

# comment on def of $nprint: $ndata[1]=nstep[1]/$nprint;
# comment on def of $nprint: $ndata[2]=nstep[2]/$nprint;
#
# Total number of trajectories to be harvested
$maxgen=300;
#
# starting value of iteration: zero for fresh run, non-zero if re-start
# after a crash
$generation=0;
#
# name of trajectory files and output files for the starting iteration
$workfile[1]="gn_000_xpr1";
$workfile[2]="gn_000_xpr2";
#
# ***************** END SECTION 1 ******************************


# ***************** SECTION 2 **********************************
#
# changing to executable directory
chdir "$execdir";

# beginning to record log files
$countfile="race_history.log";
open (count,">>$countfile");
print count "progress report","\n";
close(count);

# initializing some file names associated with charmm and some values
system "grep 'NSTEP' $workfile[1].mdout >> $workfile[1].tmp";
system "grep 'NSTEP' $workfile[2].mdout >> $workfile[2].tmp";

$dummy1=$workfile[1] . ".tmp";
$nframe[1] = `wc -l < $dummy1`;
$nframe[1] = $nframe[1]-2;
# without below line, a line is added after nframe[1].
$nframe[1]=($nframe[1])/1;
$dummy2=$workfile[2] . ".tmp";
$nframe[2] = `wc -l < $dummy2`;
$nframe[2] = $nframe[2]-2;
$nframe[2]=($nframe[2])/1;
$ndata[1]=$nframe[1];
$ndata[2]=$nframe[2];
#below line: Just to check
print "\n\n\n$ndata[1], $ndata[2]\n\n\n";

# preparing to run bead_xxx_ck to figure out the number of beads
# the number of beads is worked out according to the configurational
# bias algorithm worked out in the primary reference

$countfile="shckfls";
open (count2,">$countfile");
print count2 "2","\n";
foreach $index (1..2){
    $dummy = $workfile[$index];
    print count2 $dummy, "  ", $ndata[$index],"\n";
}
close(count2);
sleep 0.1;
system "$beadexec";
sleep 1;
$dummy= $workfile[2] . "_bdd";
open(FILE, "<$dummy");
@bead_old = <FILE>;
close(FILE);
system "rm shckfls";
sleep 1;
#Need to delete *.tmp, otherwise being accumulated.
foreach $index (1..2) {
    $dummy = $execdir . "/$workfile[$index].tmp";
    system "rm $dummy";
    sleep 1;
}

# assigning the current working files
$countfile="current_file";
open (count,">$countfile");
print count $workfile[1],"\n";
print count $workfile[2],"\n";
close(count);
# total number of frames = 2 * nsteps / nprint
$nframe[3] = $nframe[1] + $nframe[2];
# maximum_delta_shift = fraction of nframe[3] is set in section 1
$mxdlsh=int($maximum_delta_shift*$nframe[3])+1;
#
# initializing acceptance indices
$accept=1;
$naccept=0.0;
#$acceptancerate=$naccept/$generation;
$acceptancerate=0.0;
#
# Reset dis_fact to 1.0 if running BOLAS
if ($tps_or_bolas == 2) {
    $dis_fact=1.0;
}
#
# Definition: accept=1 implies accept,
# accept=0 implies reject
#
# Record sampling and log files
$countfile="sampling";
open (count,">>$countfile");
print count $workfile[1],"\n";
print count $workfile[2],"\n";
print count $nframe[1], " ",$nframe[2],"\n";
close(count);
#
$countfile="race_history.log";
open (count4,">>$countfile");
print count4 "generation ",$generation,"\n";
print count4 "acceptance rate= ", $acceptancerate,"\n";
print count4 "displacement step= ", $dis_step,"\n";
print count4 "\n";
close(count4);
sleep 0.1;
#
# Write generation_status and end files
# to terminate a run anytime type "cp end generation_status"
# at the command line and the program will exit after the current
# iteration
#
$countfile="generation_status";
open (count,">$countfile");
print count "progenesis";
close(count);
#
$countfile="end";
open (count,">$countfile");
print count "decimation";
close(count);
$decimation="decimation";
#
sleep 0.1;
#
# ************************** END SECTION 2 *************************
#
# ************************** BEGIN SECTION 3 ***********************
# begin iteration loop
#
while ($generation <= $maxgen) {

#   iteration index
    $generation+=1;
    print "generation= ",$generation,"\n";
#   irand2 and irand are variables associated with random number generator
    $irand2=`echo $irand | $random`;
    $irand=-int($irand2*100000.0);

# invoke SHOOTING/SHIFTING MOVE

    # STEP1: choose between shooting and shifting
    $irand2=`echo $irand | $random`;
    $irand=-int($irand2*100000.0);
    $dummy = int($irand2*(1+$FRQshift))+1;
    $dstp = $dis_step;
    if ($dummy > $FRQshift) {
	# then shifting chosen
	$dstp = 0.0;
	$irand2=`echo $irand | $random`;
	$irand=-int($irand2*100000.0);
	$dlsh=int($irand2*$mxdlsh)+1;
	$irand2=`echo $irand | $random`;
	$irand=-int($irand2*100000.0);
	$dummy=int($irand2*2)+1;
	$dlsh=((-1)**$dummy)*$dlsh;
	$countfile="race_history.log";
	open (count,">>$countfile");
	print count "SHIFTING MOVE INVOKED","\n";
	print count "delta shift = ",$dlsh,"\n";
	close(count);
    } else {
	# shooting chosen
	$dlsh=0;
	$countfile="race_history.log";
	open (count,">>$countfile");
	print count "SHOOTING MOVE INVOKED","\n";
	print count "delta shoot = ",$dstp,"\n";
	close(count);
    }
    # END STEP 1
    print "Step 1 done!","\n";

    # STEP2: choose a time slice at random
    $irand2=`echo $irand | $random`;
    $irand=-int($irand2*100000.0);
    $dummy=int($irand2*$bead_old[0])+1;
    $chstp=$bead_old[$dummy];

    $countfile="race_history.log";
    open (count,">>$countfile");
    print count $chstp, " chosen tstep","\n";
    close(count);
    # END STEP 2
    print "Step 2 done!","\n";

    # STEP 3: Figure out if the trajectories collectively
    #         visit basins A and B
    foreach $index (1..2) {
	$dummy= $workfile[$index] . "_out";
	$lc[$index] = `wc -l < $dummy`;
	sleep 0.1;
	$lc[$index]=$lc[$index]/(2+$nopara);
	open(FILE, "<$dummy");
	@outcheck_stuff = <FILE>;
	close(FILE);
	$reaca[$index]=0;
	$reacb[$index]=0;
	foreach $jloop (1..$lc[$index]) {
	    $dummy=(2+$nopara)*$jloop-(2+$nopara);
	    $opara = @outcheck_stuff[$dummy];
	    chomp($opara);
	    if ($opara == 1) {
		$reacb[$index]+=1;
	    } elsif ($opara == -1) {
		$reaca[$index]+=1;
	    }
	}
    }

    $connecta=0;
    $connectb=0;
    $tps_connect=1;
    $bolas_connect=1;
#   Have the old trajectories visited A?
    if ($reaca[1] > 0) {
	$connecta=1;
    } elsif ($reaca[2] > 0) {
	$connecta=1;
    }

#   Have the old trajectories visited B?
    if ($reacb[1] > 0) {
	$connectb=1;
    }elsif ($reacb[2] > 0) {
	$connectb=1;
    }

    $connect = $connecta + $connectb;
    if ($connect == 2) {
#   old trajectories have collectively visited A and B
	$tps_connect=0;
    }

    print "connect: $connecta $connectb $connect","\n";

#   Have the old trajectories visited the bead domain, i.e.,
#   the boundaries set by the bead window in bead_xxx_ck?
    if ($bead_old[0] > 0) {
#   answer is yes
	$bolas_connect=0;
    }

#   Check if the initial trajectories are valid ones for TPS/BOLAS runs
    $dummy=1;
    if ($tps_or_bolas == 1) {
	if ($tps_connect == 0) {
	    $dummy=0;
	}
    }
    if ($tps_or_bolas == 2) {
	if ($bolas_connect == 0) {
	    $dummy=0;
	}
    }

#   if dummy = 0, trajectories are valid
    if ($dummy == 0) {
	# valid trajs, pl script continues
	$countfile="race_history.log";
	open (count,">>$countfile");
	print count "Valid trajs- script proceeds","\n";
	close(count);
    } else {
	# invalid trajs, pl script exiting
	$countfile="race_history.log";
	open (count,">>$countfile");
	print count "Invalid trajs- script exiting","\n";
	close(count);
	exit(1);
    }
    # END STEP 3
    print "Step 3 done!","\n";


    # STEP 4: figure out start and end points of chain and
    # where the chosen step lies
    # for indexing, always consider traj_1 as beginner and traj_2
    # as finisher

    if ($chstp == $nframe[1]) {
	$irand2=`echo $irand | $random`;
	$irand=-int($irand2*100000.0);
	$chstp+=(2*int($irand2*2)-1);
    }

    if ($chstp < $nframe[1]) {
	$chfile = $workfile[1];
	$scale[1]=1;
	$scale[2]=-1;
#	$dummy = $nframe[1] - $chstp;
        $dummy = $nframe[1] +1 - $chstp;
	$nframe[1] = $chstp;
	$chstp = $dummy;
	$nframe[2]+=($chstp)-1;
    } else {
	$chfile = $workfile[2];
	$scale[1]=-1;
	$scale[2]=1;
	$dummy = $chstp - $nframe[1];
	$nframe[1] = $chstp;
	$chstp = $dummy;
	$nframe[2]+=((-1)*$chstp);
    }

    $nframe[1]=$nframe[1]+$dlsh;
    $nframe[2]=$nframe[2]-$dlsh;
    if ($nframe[1] < 40) {
	$nframe[1]=501;
	$nframe[2]=499;
    }
    if ($nframe[2] < 40) {
	$nframe[1]=501;
	$nframe[2]=499;
    }
    $ndata[1]=$nframe[1];
    $ndata[2]=$nframe[2];
    # end avoid problems

    # ignore top, always equal for BOLAS
    if ($tps_or_bolas == 2) {
	$nframe[1] = 500;
	$ndata[1] = 500;
	$ndata[2] = 500;
	$nframe[2] = 500;
    }
    # END STEP 4
    print "Step 4 done!","\n";
    #$chstp *500 to align with NAMD output format
    $countfile="race_history.log";
    open (count,">>$countfile");
    print count "NFRAME ", $nframe[1], " ", $nframe[2],"\n";
    print count "CHFILE ", $chfile,"\n";
    print count "CHSTP ", $chstp,"\n";
    close(count);

#   MOD(3/20) change chstp below as selected
    $selected = $chstp * $nprint;

# Amber ptraj converting binary DCD to ASCII
    $countfile="ptraj-xv" . $index . ".inp";
    open (count,">$countfile");
    print count "#################  AMBER PTRAJ  ################\n\n";
    print count "trajin  ",$chfile,".nc  ",$chstp,"  ",$chstp," \n\n";
    print count "trajout TEMP.ncrst ncrst ","\n";
    print count "\n";
    close(count);
    sleep 1;

#    $countfile="ptraj-v" . $index . ".inp";
#    open (count,">$countfile");
#    print count "#################  AMBER PTRAJ  ################\n\n";
#    print count "trajin  ",$chfile,".mdvel  ",$chstp,"  ",$chstp," \n\n";
#    print count "trajout TEMP.vel ","\n";
#    print count "\n";
#    close(count);
#    sleep 1;

#    $countfile="ptraj-xv.scr";
#    open (count,">$countfile");
#    print count "#!/bin/bash","\n\n";
#    print count "#SBATCH --job-name=",$countfile,"\n";
#    print count "#SBATCH --error=job.%J.err","\n";
#    print count "#SBATCH --output=job.%J.out","\n";
#    print count "#SBATCH --nodes=1","\n";
#    print count "#SBATCH --time=00:30:00","\n";
#    print count "module load amber/openmpi/intel/22.00";
#    print count "/share/apps/amber/22.00/openmpi/intel/bin/cpptraj ",$topology," <ptraj-x.inp> ptraj-x.outt\n\n";
#    print count "/share/apps/amber/22.00/openmpi/intel/bin/cpptraj ",$topology," <ptraj-v.inp> ptraj-v.outt\n\n";

    system "/share/apps/amber/22.00/openmpi/intel/bin/cpptraj -p $topology -i ptraj-xv.inp -o ptraj-xv.outt";
#    system "/share/apps/amber/22.00/openmpi/intel/bin/cpptraj $topology <ptraj-v.inp> ptraj-v.outt";
    $dummy = "ptraj-v.inp_done3";
    system "touch $dummy";
#    print count "touch $dummy;","\n\n";
#    close(count);
#    sleep 0.1;
#    system "chmod u+x $countfile";
#    sleep 2;
#    system "$communication $countfile";



    $filecount=0;
    while ($filecount < 1) {
      sleep 10;
      $filecount=0;
      $dummy= "ptraj-v.inp_done3";
      $filecount+=1 if -e $dummy;
    }
    print "current filecount: $filecount\n";

    system "module load netcdf-c/intel/4.7.4";
    system "/share/apps/netcdf-c/4.7.4/intel/bin/ncdump -f f TEMP.ncrst > TEMP.cdl";
    my $cmd = q{sed -n '/coordinates =/,/velocities =/p' TEMP.cdl > TEMP.crd.cdl};
    system $cmd;
    #print "extract coordinates","\n";
    my $cmd = q{sed -n '/velocities =/,/cell_spatial =/p' TEMP.cdl > TEMP.vel.cdl};
    system $cmd;
    #print "extract velocities","\n";
    system "head -n -2 TEMP.crd.cdl > crdTEMP.cdl && mv crdTEMP.cdl TEMP.crd.cdl";
    system "head -n -2 TEMP.vel.cdl > velTEMP.cdl && mv velTEMP.cdl TEMP.vel.cdl";

# Above Amber ptraj generate ascii TEMP.x TEMP.vel. Use convert code, rewrite them in pdb format.
#    $countfile="shckfls";
#    open (count5,">$countfile");
#    $dummy = $chfile;
#    print count5 $dummy,"  ",$selected,"\n";
#    close(count5);

#    system "$goconvert";
#    system "$goconvert-vel";
#    sleep 2;
#    system "rm shckfls";
#    sleep 2;

#    $filecount=0;
#    while ($filecount < 1) {
#	sleep 2;
#	$filecount=0;
#	foreach $index (1) {
#	  $dummy=$chfile . ".PERT.vel";
#	  $filecount+=1 if -e $dummy;
#	}
#    }

    # Applying velocity perturbtion #

    # Rename selected *.vel file to *.PERT.vel
#    $dummy= $chfile.".vel";
#    $dummy2 = "./".$chfile.".PERT.vel";
#    system "cp $dummy $dummy2";

    # Prepare shckfls for velocity perturbation
    LOOP: {
    $countfile="shckfls";
    open (count5,">$countfile");
    print count5 "2 ","\n";
    foreach $index (1..2){
#      $dummy = $chfile;
#      $irand2=`echo $irand | $random`;
#      $irand=-int($irand2*100000.0);
      $idummy=(-1)*$irand+202*$phase*int(($index+1)/2);
      $kdummy = $scale[$index];
#      print("$idummy \n");
      print count5 "TEMP"," ",$dstp," ",$idummy," ",$selected," ",$kdummy,"\n";
    }
    close(count5);

    #FORTRAN code for velocity perturbation
    system "$goperturb";
    sleep 2;
    system "rm shckfls";
    system "rm *.inp";
    system "rm *.outt";

    $filecount=0;
    while ($filecount < 1) {
	sleep 2;
	$filecount=0;
	foreach $index (1) {
	  $dummy="TEMP.done2.vel.cdl";
	  $filecount+=1 if -e $dummy;
	}
    }

    $countfile="raceBookKeep.log";
    open (count7,">>$countfile");
    print count7 "generation",$generation,"\n";
    #print count7 "Coordinates     ",$chfile,".PERT.coor\n";
    #print count7 "Velocities      ",$chfile,".PERTdone",$index,".vel\n";
    #print count7 "extendedSystem  ",$chfile,".xsc\n";
    close(count7);

    my $cmd = q{sed -i '1 i\ velocities =' TEMP.done1.vel.cdl};
    system $cmd;
    #print "add velocities = line","\n";

    my $cmd = q{sed -i -e '$a\ \n cell_spatial = "abc";  // cell_spatial(3)' TEMP.done1.vel.cdl};
    system $cmd;
    #print "add cell_spatial = line","\n";

    system "cp TEMP.cdl TEMP.1.cdl";

    my $cmd = q{sed -i -e '/velocities =/,/cell_spatial =/!b' -e '/cell_spatial =/!d;r TEMP.done1.vel.cdl' -e 'd' TEMP.1.cdl};
    system $cmd;
    #print "update velocities","\n";

    my $cmd = q{sed -i '1 i\ velocities =' TEMP.done2.vel.cdl};
    system $cmd;
    #print "add velocities","\n";

    my $cmd = q{sed -i -e '$a\ \n cell_spatial = "abc";  // cell_spatial(3)' TEMP.done2.vel.cdl};
    system $cmd;
    #print "add cell_spatial","\n";

    system "cp TEMP.cdl TEMP.2.cdl";

    my $cmd = q{sed -i -e '/velocities =/,/cell_spatial =/!b' -e '/cell_spatial =/!d;r TEMP.done2.vel.cdl' -e 'd' TEMP.2.cdl};
    system $cmd;
    #print "update velocities 2","\n";

    system "/share/apps/netcdf-c/4.7.4/intel/bin/ncgen -o TEMP.convert1.ncrst TEMP.1.cdl";
    system "/share/apps/netcdf-c/4.7.4/intel/bin/ncgen -o TEMP.convert2.ncrst TEMP.2.cdl";

    # STEP 5: tps/BOLAS input and script file preparation
    foreach $index (1..2) {
	$countfile="job_pr" . $index . ".inp";
#	open (count,">$countfile");
	$charmfile[$index]=$countfile;
	system "cat $tpsbarefile >> $countfile";


	system "sed -i 's/set_nprint/$nprint/g' $countfile";

	$dummy = $nframe[$index]*$nprint;
	$trajnstep[$index]=$dummy;

	system "sed -i 's/set_nsteps/$dummy/g' $countfile";

        $asfile[$index]="job_as" . $index . ".inp";
#        $asfile2[$index]="job_as2" . $index . ".inp";
	if ($generation < 10) {
	    $charmoutput[$index] = $job_file . "00" . $generation . "_xpr" . $index;
	} elsif  ($generation < 100) {
	    $charmoutput[$index] = $job_file . "0" . $generation . "_xpr" . $index;
	} else {
	    $charmoutput[$index] = $job_file . $generation . "_xpr" . $index;
	}
	# PREPARE first part of charmm input file, i.e., set variables
#	print count "#################  TPS/BOLAS with NAMD  ################\n\n";
#	print count "structure          ", $execdir,"/",$topology,"\n";
#	print count "coordinates        ", $execdir,"/oxoOPENsyn-sol.pdb\n";
#	print count "\nif {1} {\n";
#	print count "set inputname      ",$workfile[$index],".",$chstp,"\n";
#	print count "Coordinates     ",$chfile,".PERT.coor\n";
#	print count "Velocities      ",$chfile,".PERTdone",$index,".vel\n";
#	print count "extendedSystem  ",$chfile,".xsc\n";
#       print count "}\n";
#	print count "\n";
#	print count "set outputname     ",$charmoutput[$index],"\n";
#	print count "restartname        ",$charmoutput[$index],"\n";
#	print count "firsttimestep      0\n";
#	print count "#Input\n";
#	print count "paraTypeCharmm	    on\n";
#	print count "parameters          ",$execdir,"/par_all36_na.prm\n";
#       print count "parameters          ",$execdir,"/toppar_water_ions.str\n";
#	$dummy = $nframe[$index]*$nprint;
#	$trajnstep[$index]=$dummy;
#	print count "set numsteps        ",$dummy,"\n";
#       MOD(3/20) noutput is changed to nprint
#	print count "set noutput        ",$nprint;
#	print count "\n";
#	print count "\n";
#	print count "set scale = ", $scale[$index],"\n";
#	$dummy = $nframe[$index]*$nprint;
#	print count "set nsteps = ", $dummy,"\n";
#	$dummy = $chstp*$nprint;
#	print count "set chstp = ", $dummy,"\n";
#	print count "set dstp = ", $dstp,"\n";
#	print count "set timestep = ", $timestep,"\n";
#       changing seed for each generation
#	$idummy=(-1)*$irand+202*$phase*int(($index+1)/2);
#	print count "set seed = ", $idummy,"\n";
#	print count "set topfile = ", $topfile,"\n";
#       print count "set parfile = ", $parfile,"\n";
#	print count "set ndata = ", $ndata[$index],"\n";
#	print count "\n";
#	print count "\n";
#	close(count);
#	sleep 1;
	# END prepare
#        system "cat $charmfile[$index] $tpsbarefile2 > $asfile[$index]";
#	system "cat $tpsbarefile >> $charmfile[$index]";
    }

    foreach $index (1..2) {
	$countfile="job_as" . $index . ".inp";
	open (count,">$countfile");
	# PREPARE Amber PTRAJ FILE
	print count "#################  AMBER PTRAJ  ################\n\n";
	print count "trajin  ",$charmoutput[$index],".nc\n\n";
	#print count "reference step4.1_equilibration.rst7 [model1]\n";
	#print count "reference step4.1_equilibration.next.rst7 [model2]\n";
	#print count "rmsd rmsd_nonH1 ", "\@1-2465&!\@H= ","first out ",$charmoutput[$index],".dat1 mass\n";
	#print count "rmsd rmsd_nonH2 ", "\@1-2465&!\@H= ","ref [model2] out ",$charmoutput[$index],".dat2 mass\n";
	print count "distance d1 ","\:3\@P ","\:62\@P ", "out ",$charmoutput[$index],".dat1\n";
	print count "distance d2 ","\:20\@P ","\:28\@P ", "out ",$charmoutput[$index],".dat2\n";
	print count "distance d3 ","\:4\@O4' ","\:76\@N1,C2,O2,N3,C4,O4,C5,C6 ", "out ",$charmoutput[$index],".dat3\n";
	#print count "distance d01 ", "\:2\@N1 ","\:23\@N3 ", "out ", $charmoutput[$index], ".dat2\n";
	#print count "distance d11 ", "\:2\@N2 ","\:23\@O2 ", "out ", $charmoutput[$index], ".dat3\n";
	#print count "distance d21 ", "\:2\@O6 ", "\:23\@N4 ","out ", $charmoutput[$index], ".dat4\n";

#	print count "dihedral DHR2 ","\@5163 ","\@5161 ","\@5158 ","\@5155 ","out ",$charmoutput[$index],".dat5\n";
#	print count "dihedral DHR3 ","\@5393 ","\@5395 ","\@5398 ","\@5405 ","out ",$charmoutput[$index],".dat6\n";
	print count "\n";
  #print count "trajout tempTRAJ-",$index,".nc  amber \n";
	close(count);
	sleep 1;
	# END prepare
	#system "cat $tpsbarefile2 >> $asfile[$index]";
    }

#    foreach $index (1..2) {
#	$countfile="job_as2" . $index . ".inp";
#	open (count,">$countfile");
	# PREPARE Amber PTRAJ FILE
#	print count "#################  AMBER PTRAJ  ################\n\n";
#	print count "trajin  tempTRAJ-",$index,".dcd \n";
#        print count "reference ",$reference,"\n\n";
#        print count "rmsd reference out ",$charmoutput[$index],".dat8 :307-327 nofit";
#	close(count);
#	sleep 1;
	# END prepare
#        system "cat $tpsbarefile2 >> $asfile2[$index]";
#    }

    # *******************************************
    # prepare PBS script for NAMD
    foreach $index (1..2) {
	$countfile=$jobname . $index . ".scr";
        $startfile="TEMP.convert".$index.".ncrst";
	open (count,">$countfile");
	print count "#!/bin/bash","\n\n";
	print count "#SBATCH --job-name=",$countfile,"\n";
	print count "#SBATCH --error=job.%J.err","\n";
	print count "#SBATCH --output=job.%J.out","\n";
	print count "#SBATCH --nodes=1","\n";
	print count "#SBATCH --tasks-per-node=1","\n";
	#print count "#SBATCH --cores-per-socket=4","\n";
	#print count "#SBATCH --cpus-per-task=1","\n";
	print count "#SBATCH --time=2:00:00","\n";
	print count "#SBATCH --mem=10GB","\n";
  	print count "#SBATCH --gres=gpu:1","\n";
	#print count "#SBATCH --wait","\n";
	#print count "#SBATCH --requeue","\n";
	print count "module purge","\n";
	print count "module load amber/openmpi/intel/22.00","\n";
  	print count "pmemd.cuda -O -i $charmfile[$index] -p $topology -c $startfile -o $charmoutput[$index].mdout -r $charmoutput[$index].rst7 -inf $charmoutput[$index].info -x $charmoutput[$index].nc","\n";
	$dummy = $execdir . "/" . $charmfile[$index] . "_done";
	print count "wait","\n";
	print count "touch $dummy","\n";
	print count "exit","\n";
	close(count);
	sleep 0.1;
        # *******************************************
	# run the tps for this generation, i.e., two successive charrm
	# jobs for forward and backward trajectories
	system "chmod u+x $countfile";
	sleep 2;
	system "$communication $countfile >> race_history.log";
	sleep 1;
    }
    # prepare PBS script for Amber ptraj
    print "Step 5 done!","\n";

    # STEP 6
    # wait till all charmm jobs are done
    $countfile="race_history.log";
    open (count,">>$countfile");
    print count "WAITING FOR NAMD JOBS TO BE DONE","\n";
    close(count);

    $filecount=0;

############################################################
#########  KEEP FOR TEST, REMOVE '#' LATER  ################
############################################################
    while ($filecount < 2) {
	sleep 60;
	$filecount=0;
	foreach $index (1..2) {
	    $dummy=$charmfile[$index] . "_done";
	    $filecount+=1 if -e $dummy;
	}
   }
   }
    sleep 0.1;

    $traj1=$charmoutput[1].".nc";
    $traj2=$charmoutput[2].".nc";

    if (not((-e $traj1) and (-e $traj2))){
        goto LOOP;
    }

    # end wait
    # generating dat6 files
    sleep 0.1;
    foreach $index (1..2) {
#	$countfile=$jobname . $index . "ptraj.scr";
#	open (count,">$countfile");
#	print count "#!/bin/bash","\n\n";
#	print count "#SBATCH --job-name=",$countfile,"\n";
#	print count "#SBATCH --error=job.%J.err","\n";
#	print count "#SBATCH --output=job.%J.out","\n";
#	print count "#SBATCH --nodes=1","\n";
#	print count "#SBATCH --time=00:30:00","\n";
#	print count "module load amber/openmpi/intel/22.00";
#	print count "/share/apps/amber/22.00/openmpi/intel/bin/cpptraj ",$topology," <job_as",$index,".inp> ",$charmoutput[$index],".outt","\n\n";
	system "/share/apps/amber/22.00/openmpi/intel/bin/cpptraj -p $topology -i job_as$index.inp -o $charmoutput[$index].outt";
	$dummy = $execdir . "/" . $charmfile[$index] . "_done2";
#	print count "touch $dummy;","\n\n";
        system "touch $dummy";
#	close(count);
#	sleep 0.1;
#        system "chmod u+x $countfile";
#        sleep 2;
#        system "$communication $countfile >> race_history.log";
#        sleep 1;
    }

    $filecount=0;
    while ($filecount < 2) {
	sleep 10;
	$filecount=0;
	foreach $index (1..2) {
	    $dummy= $charmfile[$index] . "_done2";
	    $filecount+=1 if -e $dummy;
	}
    }

#    foreach $index (1..2) {
#	$countfile=$jobname . $index . "ptraj2.scr";
#	open (count,">$countfile");
#	print count "#!/bin/bash","\n\n";
#	print count "#SBATCH --job-name=",$countfile,"\n";
#	print count "#SBATCH --error=job.%J.err","\n";
#	print count "#SBATCH --output=job.%J.out","\n";
#	print count "#SBATCH --nodes=1","\n";
#	print count "#SBATCH --time=00:30:00","\n";
#	print count "/share/apps/amber/22.00/openmpi/intel/bin/cpptraj ",$topology," <job_as2",$index,".inp> ",$charmoutput[$index],".outt","\n\n";
#	$dummy = $execdir . "/" . $charmfile[$index] . "_done4";
#	print count "touch $dummy;","\n\n";
#	close(count);
#	sleep 0.1;
#       system "chmod u+x $countfile";
#        sleep 2;
#        system "$communication $countfile >> race_history.log";
#        sleep 1;
#    }
############################################################
#########  KEEP FOR TEST, REMOVE '#' LATER  ################
############################################################
#    $filecount=0;
#    while ($filecount < 2) {
#	sleep 10;
#	$filecount=0;
#	foreach $index (1..2) {
#	    $dummy= $charmfile[$index] . "_done4";
#	    $filecount+=1 if -e $dummy;
#	}
#    }

    # charmm jobs are done
    $countfile="race_history.log";
    open (count,">>$countfile");
    print count "NAMD JOBS ARE DONE FOR GEN = ", $generation,"\n";
    close(count);
    # END STEP 6
    print "Step 6 done!","\n";

    # STEP 7
    # program to run orderpara check
    $countfile="shckfls";
    open (count2,">$countfile");
    print count2 "2","\n";
    foreach $index (1..2){
	$dummy = $charmoutput[$index];
	print count2 $dummy, "  ", $ndata[$index],"\n";
    }
    close(count2);

    sleep 1;
    # run tps_xxx_ck
    foreach $index (1..2){
#        system "sed -i '/^#/d' $charmoutput[$index].dat1";
#        system "sed -i '/^#/d' $charmoutput[$index].dat2";
#        system "sed -i '/^#/d' $charmoutput[$index].dat3";
	system "echo $index | $shootexec >> race_history.log";
    }
    sleep 1;
    # wait till shootcheck is done
    $filecount=0;
    while ($filecount < 2) {
	sleep 10;
	$filecount=0;
	foreach $index (1..2) {
	    $dummy= $charmoutput[$index] . "_out";
	    $filecount+=1 if -e $dummy;
	}
    }
    sleep 1;
    # run bead_xxx_ck
    system "$beadexec";
    sleep 1;
    $dummy= $charmoutput[2] . "_bdd";
    open(FILE, "<$dummy");
    @bead_new = <FILE>;
    close(FILE);
    sleep 1;
    # end wait

    # opcheck donedone
    $countfile="race_history.log";
    open (count,">>$countfile");
    print count "OP CHECK DONE FOR GEN = ", $generation,"\n";
    close(count);
    # END STEP 7
    print "Step 7 done!","\n";

    # STEP 8
    # read from output the opara
    foreach $index (1..2) {
	$dummy= $charmoutput[$index] . "_out";
	$lc[$index] = `wc -l < $dummy`;
	sleep 1;
	$lc[$index]=$lc[$index]/(2+$nopara);
	open(FILE, "<$dummy");
	@outcheck_stuff = <FILE>;#$nprint=20;
	close(FILE);
	$reaca[$index]=0;
	$reacb[$index]=0;
	foreach $jloop (1..$lc[$index]) {
	    $dummy=(2+$nopara)*$jloop-(2+$nopara);
	    $opara = @outcheck_stuff[$dummy];
	    chomp($opara);
	    if ($opara == 1) {
		$reacb[$index]+=1;
	    } elsif ($opara == -1) {
		$reaca[$index]+=1;
	    }
	}
    }

    # decide to accept or reject new trajectories
    $connecta=0;
    $connectb=0;
    $tps_connect=1;
    $bolas_connect=1;
#   Have the new trajectories visited A?
    if ($reaca[1] > 0) {
	$connecta=1;
    } elsif ($reaca[2] > 0) {
	$connecta=1;
    }

#   Have the new trajectories visited B?
    if ($reacb[1] > 0) {
	$connectb=1;
    }elsif ($reacb[2] > 0) {
	$connectb=1;
    }

    $connect=$connecta+$connectb;
    if ($connect==2) {
#   new trajectories have collectively visited A and B
	$tps_connect=0;
    }

#   Have the new trajectories visited the bead domain, i.e.,
#   the boundaries set by the bead window in bead_xxx_ck?
    if ($bead_new[0] > 0) {
#   answer is yes
	$bolas_connect=0;
    }

#   Check if the new trajectories are valid ones for TPS/BOLAS runs
    $dummy=1;
    $accrej=1;
    if ($tps_or_bolas == 1) {
	if ($tps_connect == 0) {
	    $accrej=0;
	}
    }
    if ($tps_or_bolas == 2) {
	if ($bolas_connect == 0) {
	    $accrej=0;
	}
    }

    # TPS bolas acceptance probabilies
    $pacc = 1.0;
    if ($accrej == 0) {
	$pacc = (1.0*$bead_new[0])/(1.0*$bead_old[0]);
    }
    $irand2=`echo $irand | $random`;
    $irand=-int($irand2*100000.0);
    if ($accrej == 0) {
	if ($irand2 > $pacc) {
	    $accrej = 1;
	    $countfile="race_history.log";
	    open (count,">>$countfile");
	    print count "REJECTED BECAUSE OF METROPOLIS","\n";
	    print count $irand2, " GT ", $pacc,"\n";
	    close(count);
	}
    }

    if ($accrej == 0) {
	# accept new Traj
	$accept == 1;
	$naccept += 1.0;
        $workfile_old[1]=$workfile[1];
	$workfile_old[2]=$workfile[2];
        $workfile[1]=$charmoutput[1];
	$workfile[2]=$charmoutput[2];
	$countfile="current_file";
	open (count,">$countfile");
	print count $workfile[1],"\n";
	print count $workfile[2],"\n";
	close(count);

	$countfile="sampling";
	open (count,">>$countfile");
	print count $workfile[1],"\n";
	print count $workfile[2],"\n";
	print count $nframe[1], " ",$nframe[2],"\n";
	close(count);

	$countfile="sampling-acc";
	open (count,">>$countfile");
	print count $workfile[1],"\n";
	print count $workfile[2],"\n";
	print count $nframe[1], " ",$nframe[2],"\n";
	close(count);

	$countfile="race_history.log";
	open (count,">>$countfile");
	print count "TRAJS ACCEPTED FOR GEN = ", $generation,"\n";
	close(count);

	# update bead indices
	foreach $index (1..$bead_old[0]){
	    $bead_old[$index]=0;
	}
	$bead_old[0]=$bead_new[0];
	foreach $index (1..$bead_new[0]){
	    $bead_old[$index]=$bead_new[$index];
	}
	foreach $index (1..$bead_new[0]){
	    $bead_new[$index]=0;
	}
	$bead_new[0]=0;
	# end update bead indices
        #delete generation-1 dcd* files to save space.
#        if ($tps_or_bolas == 2) {
##          foreach $index (1..2) {
##              $generation2 = $generation-1;
##              if ($generation2 < 10) {
##                  $charmoutput2[$index] = $job_file . "00" . $generation2 . "_xpr" . $index;
##              } elsif  ($generation2 < 100) {
##                  $charmoutput2[$index] = $job_file . "0" . $generation2 . "_xpr" . $index;
##              } else {
##                  $charmoutput2[$index] = $job_file . $generation2 . "_xpr" . $index;
##              }
##           }
#           $dummy1 = $execdir . "/" . $workfile_old[1] . ".dcd*";
#           $dummy2 = $execdir . "/" . $workfile_old[2] . ".dcd*";
#           system "rm $dummy1";
#           sleep 1;
#           system "rm $dummy2";
#           sleep 1;
#        }
    } else {
	#reject new Traj
	$accept == 0;
	#$ workfiles = workfiles old, no updates necessary
	$countfile="current_file";
	open (count,">$countfile");
	print count $workfile[1],"\n";
	print count $workfile[2],"\n";
	close(count);
	system "grep 'NSTEP' $workfile[1].mdout >> $workfile[1].tmp";
	system "grep 'NSTEP' $workfile[2].mdout >> $workfile[2].tmp";
	$dummy1=$workfile[1] . ".tmp";
	$nframe[1] = `wc -l < $dummy1`;
	# without below line, a line is added after nframe[1].
        $nframe[1] = $nframe[1]-2;
	$nframe[1]=($nframe[1])/1;
	$dummy2=$workfile[2] . ".tmp";
	$nframe[2] = `wc -l < $dummy2`;
	$nframe[2] = $nframe[2]-2;
        $nframe[2]=($nframe[2])/1;
#	$dummy=$workfile[1] . ".ene";
#	$nframe[1] = `wc -l < $dummy`;
#	$nframe[1]=($nframe[1]-5)/4;
#	$dummy=$workfile[2] . ".ene";
#	$nframe[2] = `wc -l < $dummy`;
#	$nframe[2]=($nframe[2]-5)/4;
	$nframe[3] = $nframe[1] + $nframe[2];
	$countfile="sampling";
	open (count,">>$countfile");
	print count $workfile[1],"\n";
	print count $workfile[2],"\n";
	print count $nframe[1], " ",$nframe[2],"\n";
	close(count);

	$countfile="race_history.log";
	open (count,">>$countfile");
	print count "TRAJS REJECTED FOR GEN = ", $generation,"\n";
	close(count);
############################################################
#########  KEEP FOR TEST, REMOVE '#' LATER (DOWN) ##########
############################################################
	$dummy = $execdir . "/$workfile[1]" . $index . ".tmp";
	system "rm $dummy";
	sleep 1;
	$dummy = $execdir . "/$workfile[2]" . $index . ".tmp";
	system "rm $dummy";
	sleep 1;
        $dummy1 = $execdir . "/" . $charmoutput[1] . ".nc";
        $dummy2 = $execdir . "/" . $charmoutput[2] . ".nc";
	system "rm $dummy1";
        sleep 1;
	system "rm $dummy2";
        sleep 1;
    }
#remove unwanted files
    foreach $index (1..2) {
	$dummy = $execdir . "/job_pr" . $index . ".inp_done";
	system "rm $dummy";
	sleep 1;
	$dummy = $execdir . "/job_pr" . $index . ".inp_done2";
	system "rm $dummy";
	sleep 1;
#	$dummy = $execdir . "/job_pr" . $index . ".inp_done4";
#	system "rm $dummy";
#	sleep 1;
	$dummy = $execdir . "/job_pr" . $index . ".inp";
	system "rm $dummy";
	sleep 1;
	$dummy = $execdir . "/job_as" . $index . ".inp";
	system "rm $dummy";
	sleep 1;
    }
#    $dummy=$jobname . "?.scr";
#    system "rm $dummy";
#    $dummy=$jobname . "?ptraj.scr";
#    system "rm $dummy";
    system "rm shckfls";
    system "rm *_done3";
#    system "rm slurm-*.out";
    system "rm *.outt";
    system "rm TEMP* ";
#    system "rm job.*";

############################################################
#########  KEEP FOR TEST, REMOVE '#' LATER (UP) ############
############################################################
    $acceptancerate=$naccept/$generation;

# modify $dis_step to target acceptance rate (only for TPS)
    $dummy = $dis_step;
    if ($acceptancerate <= $acc_desire) {
	$dummy = $dis_step/$dis_fact;
    }
    if ($acceptancerate > $acc_desire) {
	$dummy = $dis_step*$dis_fact;
    }
    $dis_step = $dummy;
# end modify dis_step

    $countfile="race_history.log";
    open (count4,">>$countfile");
    print count4 "generation",$generation,"\n";
    print count4 "acceptance rate=", $acceptancerate,"\n";
    print count4 "displacement step=", $dis_step,"\n";
    print count4 "\n";
    close(count4);
    sleep 2;

    # update number of processors
    $countfile="nproc";
    open (count,$countfile);
    $counter=<count>;
    close(count);
    chomp($counter);
    $nproc=$counter;

    # CHECK TO SEE THE NEED FOR TERMINATION
    $countfile="generation_status";
    open (count,$countfile);
    $counter=<count>;
    close(count);
    chomp($counter);
    if ($counter eq $decimation) {
	$countfile="race_history.log";
	open (count4,">>$countfile");
	print count4 "DECIMATION INVOKED","\n";
	close(count4);
	exit(1);
    }
    print "Step 8 done!","\n";

}
