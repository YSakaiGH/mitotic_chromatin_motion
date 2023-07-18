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

set para "_loop_ext=1_cond_del=0_box_const=0"


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
set after_loop_time  10
set after_cond_del_time 10
set msd_calculation_time 100

puts "after_loop_time(sec) = [expr $time_step*$n_step*$after_loop_time/1000]"
puts "after_cond_del_time(sec) = [expr $time_step*$n_step*$after_cond_del_time/1000]"
puts "msd_calculation_time(sec) = [expr $time_step*$n_step*$msd_calculation_time/10000]"



#####################
#  Particles Setup  #
#####################
set l_poly 5000  ;  set l_loop 50 ; set n_loop [expr $l_poly/$l_loop]
set n_cond $n_loop 
set n_poly 1
set n_part0 [expr $l_poly + $n_cond] ; set n_part [expr $n_part0*$n_poly]

# Simulation Box
set rad_sphere 23.81
set l_box [expr 2.0*$rad_sphere]


###########################
#  Particle Interactions  #
###########################

########################
#  WCA ; repulsive LJ  #
########################
inter 0 0 lennard-jones 1.0 1.0 1.12246 0.25 0.0
inter 0 1 lennard-jones 1.0 1.0 1.12246 0.25 0.0
inter 1 1 lennard-jones 1.0 1.0 1.12246 0.25 0.0


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
inter 1 10 lennard-jones 1.0 1.0 1.12246 0.25 0.0
inter 2 10 lennard-jones 1.0 1.0 1.12246 0.25 0.0

constraint cylinder center 0.0 0.0 0.0 axis 0.0 0.0 1.0 radius 16.0 length 40.0 direction -1 type 10
}



#######################
#  MSD  Calculations  #
#######################
proc msd {para l_poly ti} { 
    set msd_file [open "MSD0$para" "a"]
    set msd [analyze g123 0 1 $l_poly]  ;  puts $msd_file "$ti  [lindex $msd 0]*0.0004"  ;  close $msd_file  }


##########################
#   MSD of peri & core   #
##########################
proc msd2 {para l_poly ti} { 
    set cfile [open "cfile$para"]
    gets $cfile conf

    set pid 0  ;  set msd0 0.0  ;  set msd1 0.0  ;  set msd2 0.0
    while { $pid < $l_poly } { gets $cfile conf
	set posi [ part $pid print pos ]

#  periphery : axis = 1 : 1
	if { fmod($pid,50) <= 12 || fmod($pid,50) >= 38 } { set msd0 [ expr [veclen [vecsub $posi $conf]]**2 + $msd0 ] }
	if { fmod($pid,50) >= 13 && fmod($pid,50) <= 37 } { set msd1 [ expr [veclen [vecsub $posi $conf]]**2 + $msd1 ] }
        incr pid }

    set msd0 [expr $msd0 / 2500]  ;  set msd1 [expr $msd1 / 2500]
    set msd_file [open "MSD1$para" "a"]  ;  puts $msd_file "$ti/100.0  $msd0*0.0004"  ;  close $msd_file
    set msd_file [open "MSD2$para" "a"]  ;  puts $msd_file "$ti/100.0  $msd1*0.0004"  ;  close $msd_file
}




#########################
#   Starting Condition  #
#########################
set conffile [open "polymer.init"]
gets $conffile conf

######################
#  chromatin chains  #
######################
set pid 0
while { $pid < $l_poly } { gets $conffile conf
set x1 [lindex $conf 0]  ;  set y1 [lindex $conf 1]  ;  set z1 [lindex $conf 2]
part $pid pos [expr 0.8*$x1] [expr 0.8*$y1] [expr 0.8*$z1] type 0
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


set vtf_file [open "test.vtf" "w"]  ;  writevsf $vtf_file  ;  close $vtf_file


###########################
##   Condensin Adhesion  ##
###########################
set pid 0
while { $pid < $n_cond } { set pid_c [ expr $pid + $l_poly ]
    part [expr $pid + $l_poly] bond $HARM_c [expr int(($pid+0.5)*$l_loop)]
    if { $n_poly==2 } { part [expr $n_part0 + $pid + $l_poly] bond $HARM_c [expr int( $n_part0 + ($pid+0.5)*$l_loop )] }
    incr pid }


############
#  Warmup  #
############
puts "Warmup start"
set min 0  ;  set cap 10  ;  set rcap 1000
while { $cap <= $rcap } { inter forcecap $cap  ;  integrate 10  ;  set min [analyze mindist]  ;  incr cap 10 }

set vtf_data [open "test000.dat" "w"] ; writevcf $vtf_data ; close $vtf_data



#########################
##    Loop Extrusion   ##
#########################
if { $loop_extrusion == 1 } {

set loop_step 1000

set loop_i 1
while { $loop_i <= [expr $l_loop/2 - 1] } {


#############################
#   condensin step forward  #
#############################
    set pid 0
    while { $pid < $n_cond } { 
	set half [expr int( ($pid + 0.5)*$l_loop )]
	set pls [expr $half + $loop_i]  ;  set mns [expr $half - $loop_i]
	part [expr $pid + $l_poly] bond delete
	part [expr $pid + $l_poly] bond $HARM_c $mns  ;  part [expr $pid + $l_poly] bond $HARM_c $pls
	incr pid }
#############################

    integrate $loop_step

    set vtf_file [open "test.vtf" "a"] ; writevcf $vtf_file ; close $vtf_file
    puts "loop_length = $loop_i"
    
    incr loop_i }

######################################
#  condensin base position of loops  #
######################################
set pid 0
while { $pid < [expr $n_cond - 1] } {
	part [expr $pid + $l_poly] bond delete
	part [expr $pid + $l_poly] bond $HARM_c [expr $pid*$l_loop]
	part [expr $pid + $l_poly] bond $HARM_c [expr ($pid+1)*$l_loop]
	incr pid  }

integrate $loop_step


set vtf_file [open "test.vtf" "a"] ; writevcf $vtf_file ; close $vtf_file
puts "loop_length = $loop_i"

}



####################
#  Thermarization  #
####################
set i 1
while { $i <= $after_loop_time } { integrate $n_step  ;  puts "after_cond_del_time = $i"
    set vtf_file [open "test.vtf" "a"] ; writevcf $vtf_file ; close $vtf_file
    incr i  }


######################
#  Condensin Delete  #
######################
if { $condensin_deletion == 1 } {
set pid 0
while { $pid < $n_cond } { part [expr $pid + $l_poly] delete  ;  incr pid }
}


####################
#  Thermarization  #
####################
set i 1
while { $i <= $after_cond_del_time } { integrate $n_step  ;  puts "after_cond_del_time = $i"
    set vtf_file [open "test.vtf" "a"] ; writevcf $vtf_file ; close $vtf_file
    incr i  }


####################
#  Thermalization  #
####################
analyze g123 -init 0 1 $l_poly
set data [open "cfile$para" "w"] ; writevcf $data ; close $data
msd $para $l_poly 0
msd2 $para $l_poly 0

set i 1
puts "MSD calculation"
while { $i <= $msd_calculation_time } { integrate [expr $n_step/10]  ;  puts "MSD calculation time = $i"
    set vtf_file [open "test.vtf" "a"] ; writevcf $vtf_file ; close $vtf_file
    set vtf_data [open "test$i.dat" "w"] ; writevcf $vtf_data ; close $vtf_data
    msd $para $l_poly $i
    msd2 $para $l_poly $i
    incr i }

################################
#  Store final configurations  #
################################
set vtf_final [open "final$para" "w"]
writevcf $vtf_final
close $vtf_final
