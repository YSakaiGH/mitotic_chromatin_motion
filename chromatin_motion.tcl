###########################################################
#            Simulation for Mitotic chromosome            #
#                     by using ESPReSo                    #
###########################################################


###########################
#  Simulation Conditions  #
###########################
set loop_extrusion 1
set condensin_deletion 0
set box_constraint 0

#set para "_loop1_condel0_const0"
set para ""


##############
#  MD Setup  #
##############
# Simulation setup
set tempe 1.0  ;  set gamma 1.0
set time_step 0.01  ;  set skin 5
set n_step  10000

setmd box_l  50.0 50.0 50.0
setmd periodic  0 0 0
setmd time_step $time_step
setmd skin $skin

integrate set nvt
thermostat langevin $tempe $gamma


####################
#  Time Evolution  #
####################
set after_loop_time  100
set after_cond_del_time 100
set msd_calculation_time 100

puts "######################################################################"
puts "after_loop_time(sec) = [expr $time_step*$n_step*$after_loop_time/1000]"
puts "after_cond_del_time(sec) = [expr $time_step*$n_step*$after_cond_del_time/1000]"
puts "msd_calculation_time(sec) = [expr $time_step*$n_step*$msd_calculation_time/10000]"
puts "######################################################################"


#####################
#  Particles Setup  #
#####################
set l_poly 5000  ;  set l_loop 50 ; set n_loop [expr $l_poly/$l_loop]
set n_cond $n_loop 
set n_poly 1
set n_part0 [expr $l_poly + $n_cond] ; set n_part [expr $n_part0*$n_poly]


###########################
#  Particle Interactions  #
###########################

########################
#  WCA ; repulsive LJ  #
########################
inter 0 0 lennard-jones 1.0 1.0 1.12246 0.25 0.0

##########
#  Bond  #
##########
set HARM 10
set k_harm 1000.0  ;  set r_harm 1.0  ;  set harm [expr int($r_harm*10)]
inter $HARM harmonic $k_harm $r_harm

set HARM_c 20
set k_harmc 1000.0
inter $HARM_c harmonic $k_harmc 0.0


#####################
#  Bulk Constraint  #
#####################
if { $box_constraint == 1 } {
inter 0 10 lennard-jones 1.0 1.0 1.12246 0.25 0.0
constraint cylinder center 0.0 0.0 0.0 axis 0.0 0.0 1.0 radius 16.0 length 40.0 direction -1 type 10
}


##########################
#   Definition of MSD    #
#   mean & peri & core   #
##########################
proc msd {para l_poly tim} { 
    set cfile [open "cfile$para"]
    gets $cfile conf

    set msd0 0.0  ;  set msd1 0.0  ;  set msd2 0.0

    set pid 0
    while { $pid < $l_poly } { gets $cfile conf
	set posi [ part $pid print pos ]

	#  mean
	set msd0 [ expr [veclen [vecsub $posi $conf]]**2 + $msd0 ]

	#  periphery : axis = 1 : 1
	if { fmod($pid,50) <= 12 || fmod($pid,50) >= 38 } { set msd1 [ expr [veclen [vecsub $posi $conf]]**2 + $msd1 ] }
	if { fmod($pid,50) >= 13 && fmod($pid,50) <= 37 } { set msd2 [ expr [veclen [vecsub $posi $conf]]**2 + $msd2 ] }
        incr pid }

    set msd0 [expr $msd0 / 5000]  ;  set msd1 [expr $msd1 / 2500]  ;  set msd2 [expr $msd2 / 2500]

    set msd_file [open "MSD0$para" "a"]  ;  puts $msd_file "[expr $tim/100.0]  [expr $msd0*0.0004]"  ;  close $msd_file
    set msd_file [open "MSD1$para" "a"]  ;  puts $msd_file "[expr $tim/100.0]  [expr $msd1*0.0004]"  ;  close $msd_file
    set msd_file [open "MSD2$para" "a"]  ;  puts $msd_file "[expr $tim/100.0]  [expr $msd2*0.0004]"  ;  close $msd_file
}


#########################
#   Starting Condition  #
#########################
set conffile [open "polymer.init"]


######################
#  chromatin chains  #
######################
set pid 0
while { $pid < $l_poly } { gets $conffile conf
    set x1 [lindex $conf 0]  ;  set y1 [lindex $conf 1]  ;  set z1 [lindex $conf 2]
    part $pid pos [expr $x1] [expr $y1] [expr $z1] type 0
    if { $pid > 0 } { part $pid bond $HARM [expr $pid - 1] }
    incr pid }


################
#  Condensins  #
################
set pid 0
while { $pid < $n_cond } { set posi [ part [expr $pid*$l_loop] print pos ]
    set x [lindex $posi 0] ; set y [lindex $posi 1] ; set z [lindex $posi 2]
    part [expr $pid + $l_poly] pos [expr $x+0.5] [expr $y+0.5] [expr $z+0.5] type 1
    incr pid }


set vtf_file [open "poly$para.vtf" "w"]  ;  writevsf $vtf_file  ;  close $vtf_file
set vtf_file [open "poly$para.vtf" "a"]  ;  writevcf $vtf_file  ;  close $vtf_file


###########################
##   Condensin Adhesion  ##
###########################
set pid 0
while { $pid < $n_cond } { set pid_c [ expr $pid + $l_poly ]
    part [expr $pid + $l_poly] bond $HARM_c [expr int(($pid+0.5)*$l_loop)]
    incr pid }


############
#  Warmup  #
############
puts "Warmup start"
set min 0  ;  set cap 100  ;  set rcap 1000
while { $cap <= $rcap } { inter forcecap $cap  ;  integrate 100  ;  set min [analyze mindist]  ;  incr cap 10 }
puts "Warmup finish"


#########################
##    Loop Extrusion   ##
#########################
if { $loop_extrusion == 1 } {

    set loop_step [expr $n_step/10]

puts "######################################################################"
set loop_len 1
while { $loop_len <= [expr $l_loop/2 - 1] } {

    set pjd 0
    while { $pjd < $n_poly } {

    set pid 0
    while { $pid < $n_cond } { 
	set half [expr int( ($pid + 0.5)*$l_loop + $pjd*$n_part0 )]
	set pls [expr $half + $loop_len]  ;  set mns [expr $half - $loop_len]
	part [expr $pid + $l_poly + $pjd*$n_part0] bond delete
	part [expr $pid + $l_poly + $pjd*$n_part0] bond $HARM_c $mns  ;  part [expr $pid + $l_poly + $pjd*$n_part0] bond $HARM_c $pls
	incr pid }
    incr pjd }
#########################
    integrate $loop_step

    set vtf_file [open "poly$para.vtf" "a"] ; writevcf $vtf_file ; close $vtf_file
    puts "loop_len = $loop_len"
    
    incr loop_len }

##########################
#  condensin binds edges #
##########################
set pjd 0
while { $pjd < $n_poly } {
    set pid 0
    while { $pid < [expr $n_cond - 1] } { set pid_c [ expr $pid + $l_poly + $pjd*$n_part0 ]
	part $pid_c bond delete
	part $pid_c bond $HARM_c [expr $pid*$l_loop + $pjd*$n_part0]
	part $pid_c bond $HARM_c [expr ($pid+1)*$l_loop + $pjd*$n_part0]
	incr pid  }
    incr pjd }

integrate $loop_step


set vtf_file [open "poly$para.vtf" "a"] ; writevcf $vtf_file ; close $vtf_file
puts "loop_len = $loop_len"
}



########################
#  Time after looping  #
########################
puts "######################################################################"
set i 1
while { $i <= $after_loop_time } { integrate $n_step
    puts "after_loop_time(sec) = [expr $time_step*$n_step*$i/1000]"
    set vtf_file [open "poly$para.vtf" "a"] ; writevcf $vtf_file ; close $vtf_file
    incr i  }


######################
#  Condensin Delete  #
######################
if { $condensin_deletion == 1 } {
    set pid 0
    while { $pid < $n_cond } { part [expr $pid + $l_poly] delete  ;  incr pid }
}


###################################
#  Time after condensin delation  #
###################################
puts "######################################################################"
set i 1
while { $i <= $after_cond_del_time } { integrate $n_step
    puts "after_cond_del_time(sec) = [expr $time_step*$n_step*$i/1000]"
    set vtf_file [open "poly$para.vtf" "a"] ; writevcf $vtf_file ; close $vtf_file
    incr i  }


######################
#  measure MSD step  #
######################
set data [open "cfile$para" "w"] ; writevcf $data ; close $data
msd $para $l_poly 0

puts "######################################################################"
set i 1
puts "MSD calculation"
while { $i <= $msd_calculation_time } { integrate [expr $n_step/10]
    puts "MSD calculation time(sec) = [expr $time_step*$n_step*$i/10000]"
    set vtf_file [open "poly$para.vtf" "a"] ; writevcf $vtf_file ; close $vtf_file
    msd $para $l_poly $i
    incr i }


################################
#  Store final configurations  #
################################
puts "######################################################################"
set vtf_final [open "poly$para.final" "w"]
writevcf $vtf_final
close $vtf_final
