#!/usr/bin/perl -w
open A,"tgmg.in.cc";
open B,">tgmg.cc";

# Constants needed for the mul_a routine
@c=("-m-1","-m","1-m","-1","","1","m-1","m","m+1");
@cpx=("f[-1]","f[-m]","f[1-2*m]","f[m-1]","","f[1-m]","f[2*m-1]","f[m]","f[1]");
@cpy=("f[mn-m-1]","f[mn-m]","f[mn-m+1]","f[-1]","","f[1]","f[m-1-mn]","f[m-mn]","f[m+1-mn]");
@cpxy=("f[mn-1]","","z[mn-m]","","","","z[m-1]","","(*z)");
$nl="\n                   ";

# Constants needed for the RAT computation
@elx=("l","c","r");
@ely=("d","c","u");
$no_d="} else {dl=M(0.0);dc=M(0.0);dr=M(0.0);}\n";
$no_u="} else {ul=M(0.0);uc=M(0.0);ur=M(0.0);}\n";

@fr=(0.5,1,0.5);
@descx=("left","center","right");
@descy=("Bottom","Middle","Top");

# Start of RAT computation function
$rat_start=<<EOF;
    int ci,cij;
    M dl,dc,dr,cl,cc,cr,ul,uc,ur;
EOF

# End of RAT computation function
$rat_end=<<EOF;

    // Store reciprocal of central element and update pointer
    tp[9]=1./tp[4];tp+=10;
}
EOF

# Double restriction routine
$restrict_double=<<EOF;
    V ci,ci2;
    ci=0.5*res(0,ij+m);
    cp[sm]===ci;*cp===res(0,ij)+ci;
    int i=1;ij++;
    while(i<m-1) {
        ci=0.5*res(i,ij);ci2=0.25*res(i,ij+m);
        cp[sm]+=ci2;*(cp++)+=ci+ci2;
        i++;ij++;
        ci2+=0.5*res(i,ij+m);cp[sm]===ci2;
        *cp===res(i,ij)+ci+ci2;
        i++;ij++;
    }
    if(i==m-1) {
        if(x_prd) {
            ci=0.5*res(i,ij);ci2=0.25*res(i,ij+m);
            cp[sm]+=ci2;*cp+=ci+ci2;
            cp[1]+=ci2;cp[1-sm]+=ci+ci2;
        } else {
            cp++;ci=0.5*res(i,ij+m);
            *cp===res(i,ij)+ci;cp[sm]===ci;
        }
    }
}
EOF

# Single restriction routine
$restrict_single=<<EOF;
    V ci=res(0,ij);
    *cp===ci;
    int i=1;ij++;
    while(i<m-1) {
        ci=0.5*res(i,ij);*cp+=ci;cp++;i++;ij++;
        ci+=res(i,ij);*cp===ci;i++;ij++;
    }
    if(i==m-1) {
        if(x_prd) {ci=0.5*res(i,ij);*cp+=ci;cp[1-sm]+=ci;}
        else {cp++;ci=res(i,ij);*cp===ci;}
    }
}
EOF

# Loop over the input C++ file, searching for functions that need to be
# expanded
while(<A>) {
    if(/tsub_contrib\((\d),(\d)\)/) {
        $a=$1;$b=$2;
        $r=$xp=$yp=$xyp="";
        foreach $j (0..2) {
            foreach $i (0..2) {
                $k=$i+$j*3;
                next if $k==4;
                $es=($k==0?"+*e":"+e[$k]")."*";

                if($a==$i&&$a!=1) {
                    if($b==$j&&$b!=1) {
                        $xyp.=$es.$cpxy[$k];
                    } else {
                        $xp.=$es.$cpx[$k];
                    }
                } else {
                    if($b==$j&&$b!=1) {
                        $yp.=$es.$cpy[$k];
                    } else {
                        $r.=$es."f[$c[$k]]";
                    }
                }
            }
        }
        $r=~s/^\+//g;
        $xp=~s/^\+//g;
        $yp=~s/^\+//g;
        $xyp=~s/^\+//g;
        $str=$xp ne ""?($yp ne ""?"$r+(x_prd?$xp$nl+(y_prd?$xyp:V(0.0)):(y_prd?$yp:V(0.0)))":"$r$nl+(x_prd?$xp:V(0.0))")
                  :($yp ne ""?"$r$nl+(y_prd?$yp:V(0.0))":$r);
        $str=~s/\+e\[6/$nl+e[6/ if $a==1 && $b==1;
        s/tsub_contrib\(\d,\d\)/$str/;
    } elsif (/^void tgmg_base<S,V,M>::compute_rat_bdry/) {
        rat(1);
    } elsif (/^void tgmg_base<S,V,M>::compute_rat/) {
        rat(0);
    } elsif (/^void tgmg_base<S,V,M>::r_double_line_(...)/) {
        $set=$1 eq "set";
        s/;$/ {/;
        print B;

        $_=$restrict_double;
        $set?s/===/=/g:s/===/+=/g;
    } elsif (/^void tgmg_base<S,V,M>::r_line_(...)/) {
        $set=$1 eq "set";
        s/;$/ {/;
        print B;

        $_=$restrict_single;
        $set?s/===/=/g:s/===/+=/g;
    } elsif (/^void tgmg_base<S,V,M>::r_periodic_line_add/) {
        s/;$/ {/;
        print B;
        print B "    const int sd=sm*(1-sn);\n";

        $_=$restrict_single;
        s/0.5\*res\(i,ij\)/0.25*res(i,ij)/g;
        s/=res\((.),ij\)/=0.5*res($1,ij)/g;
        s/===/+=/g;
        s/\*cp\+=ci/*cp+=ci;cp[sd]+=ci/g;
        s/cp\[1-sm\]\+=ci/cp[1-sm]+=ci;*c+=ci/g;
    }
    print B;
}

# Function to generate the RAT computation routines. If zero is passed as an
# argument, it generates the routine for computation in the bulk. If one is
# passed, it generates the routine for computation taking into account
# boundaries and periodicity.
sub rat {
    $pass=$_[0];

    # Create lists to give the displacements to nearby grid points,
    # incorporating the periodicity check on the boundary pass
    @disx=(($pass?"-(w&1?1-m:1)":"-1"),"","+1");
    @disy=(($pass?"-(w&16?m-mn:m)":"-m"),"","+m");

    # Print function header
    s/;$/ {/;
    print B;

    # Print computation header
    print B $rat_start;

#   if($pass) {
#       print B <<EOF;
#   char s[5]="HDRI";
#   putchar(s[w&3]);
#   putchar(s[(w>>2)&3]);
#   putchar(',');
#   putchar(s[(w>>4)&3]);
#   putchar(s[(w>>6)&3]);
#   putchar(' ');
#EOF
#   } else {
#       print B "    fputs(\"..,.. \",stdout);\n";
#   }
    # Set the markers used in determining whether or not an element has
    # been written to
    $cha[$_]=0 foreach 0..8;
    $sta[$_]="" foreach 0..8;

    foreach $ty (0,1,2) {foreach $tx (0,1,2) {
        $txy=$tx+3*$ty;
        $pre=$fr[$tx]*$fr[$ty];

        print B "\n    // $descy[$ty] $descx[$tx] contribution\n";

        # On boundary pass, deal with "if" statement to completely
        # exclude a contribution
        if($pass) {
            $mask=($tx==0?2:($tx==2?8:0))|($ty==0?32:($ty==2?128:0));

            if($mask) {
                print B "    if(".($tx!=1&&$ty!=1?"(w&$mask)==$mask":"w&$mask").") {\n";
                $t="        ";
            } else {
                $t="    ";
            }
        } else {
            $t="    ";
        }

        # Create indices to fine grid points if needed
        if($txy!=4) {
            $cij="cij";
            print B "${t}cij=ij$disy[$ty]$disx[$tx];";
            if($tx==1) {
                print B "\n";
                $ci="i";
            } else {
                print B "ci=i$disx[$tx];\n";
                $ci="ci";
            }
        } else {
            $ci="i";
            $cij="ij";
        }

        foreach $dj (0..2) {
            if($pass) {
                if($ty==1&&$dj==0) {
                    print B "${t}if(w&32) {\n";
                    xline(1);
                    print B "${t}} else if(w&16) {\n";
                    $pre*=2;xline(1);$pre*=0.5;
                    print B "${t}$no_d";
                } elsif($ty==1&&$dj==2) {
                    print B "${t}if(w&128) {\n";
                    xline(1);
                    print B "${t}} else if(w&64) {\n";
                    $pre*=2;xline(1);$pre*=0.5;
                    print B "${t}$no_u";
                } else {xline(0);}
            } else {

                # For bulk routine evaluate all nine matrix
                # elements, as they will be required
                xline(0);
                next;
            }
        }

        # Output code to merge common terms
        $en=$pass?"}\n":"\n";
        if($ty==0) {
            print B $t."ul+=cl;cl+=dl;uc+=cc;cc+=dc;ur+=cr;cr+=dr;\n";
        } elsif($ty==2) {
            print B $t."dl+=cl;cl+=ul;dc+=cc;cc+=uc;dr+=cr;cr+=ur;\n";
        } else {
            print B $t.($pass?"if(w&32) {":"")."cl+=dl;cc+=dc;cr+=dr;$en";
            print B $t.($pass?"if(w&128) {":"")."cl+=ul;cc+=uc;cr+=ur;$en";
        }
        if($tx==0) {
            print B $t.($ty!=0?"dr+=dc;dc+=dl;":"")
                   ."cr+=cc;cc+=cl;".($ty!=2?"ur+=uc;uc+=ul;\n":"\n");
        } elsif($tx==2) {
            print B $t.($ty!=0?"dl+=dc;dc+=dr;":"")
                   ."cl+=cc;cc+=cr;".($ty!=2?"ul+=uc;uc+=ur;\n":"\n");
        } else {
            print B $t.($pass?"if(w&2) {":"").($ty!=0?"dc+=dl;":"")."cc+=cl;".($ty!=2?"uc+=ul;$en":$en);
            print B $t.($pass?"if(w&8) {":"").($ty!=0?"dc+=dr;":"")."cc+=cr;".($ty!=2?"uc+=ur;$en":$en);
        }

        $dij=0;print B $t;
        foreach $dj (-1..1) {
            $lj=2-$ty+$dj;
            next if $lj<0 || $lj>2;
            foreach $di (-1..1) {
                $li=2-$tx+$di;
                next if $li<0 || $li>2;
                $dij=4+$di+3*$dj;

                print B ($dij==0?"*tp":"tp[$dij]")."$sta[$dij]=$ely[$lj]$elx[$li];";
                if($sta[$dij] eq "") {
                    $sta[$dij]="+";
                    $cha[$dij]=1;
                }
            }
        }
        print B "\n";

        # If needed, add bracket for "if" statement
        # excluding this term
        if($t eq "    ") {
            $cha[$_]=0 foreach 0..8;
        } else {
            $chat=0;
            $chat+=$cha[$_] foreach 0..8;
            if($chat==0) {
                print B "    }\n";
            } else {
                print B "    } else ";
                print B "{" unless $chat==1;
                foreach (0..8) {
                    print B ($_?"tp[$_]":"*tp")."=M(0.0);" if $cha[$_];
                    $cha[$_]=0;
                }
                print B $chat==1?"\n":"}\n";
            }
        }
    }}
    $_=$rat_end;
}

sub xline {
    $t.="    " if $_[0];
    print B "$t";
    if(!$pass||$tx!=1) {
        foreach $di (0..2) {
            grab($di,$dj,0);
        }
    } else {
        grab(0,$dj,1);
        print B "\n$t";grab(1,$dj,0);
        print B "\n$t";grab(2,$dj,1);
    }
    print B "\n";
    $t=substr $t,4 if $_[0];
}

sub grab {
    $fac=$pre*(($tx+$_[0])%2==0?1:0.5)*(($ty+$_[1])%2==0?1:0.5);
    $in="$ely[$_[1]]$elx[$_[0]]";
    $sfac=$fac==1?"":"$fac*";
    $ter="q.a_${in}($ci,$cij)";

    if(!$pass||$tx!=1) {
        print B "$in=$sfac$ter;";
        return;
    }

    $dfac=$fac*2;$sdfac=$dfac==1?"":"$dfac*";
    $o=$_[0]==0?"w&2?$sfac$ter:(w&1?$sdfac$ter:M(0.0))":
      ($_[0]==2?"w&8?$sfac$ter:(w&4?$sdfac$ter:M(0.0))":"$sfac$ter");
    print B "$in=$o;";
}
