
!***************************************************************************************************
!***************************************************************************************************
!   SAFT Module to:
!       1. Calculate volume
!***************************************************************************************************
!
!***************************************************************************************************
module Vol_mod
!***************************************************************************************************
!Modules
!=======    
    use Types           ! Definitions of types and double precision
    use Global          ! Important global parameters
    use Setup           ! To setup the system
    use Pressure        ! Calculate pressure
    use Ideal           ! Calculate ideal A
    use Mono            ! Calculate mono A
    use Chain           ! Calculate chain A
    use Assoc           ! Calculate assoc A
    use Ion             ! Calculate ion A
    use ChemPot         ! Calculate chem pot
!***************************************************************************************************
    implicit none
!***************************************************************************************************
    contains
!***************************************************************************************************
    function Vol_dens_g(limited) result(v_out)   
        implicit none
        
        logical, optional   ::  limited
        real(kind=DP)       ::  v_out
        real(kind=DP)       ::  lv_array(1:10), lg_array(1:10)
        real(kind=DP)       ::  v_temp1, v_temp2, p_temp1, p_temp2, d_v
        real(kind=DP)       ::  lv1, lv2, lp1, lp2, lr, g_local, grad_local
        integer             ::  i_rho, nroots, piter, froot, imusum
        integer, parameter          :: plimit1=100, plimit2=500, i_rho_limit=5000
        real(kind=DP), parameter    ::  P_crit=1.0e-6_DP,P_crit2=1.0e-4

        !set limit if required
        if(present(limited)) then
            limited = .true.
        end if
        
        nroots = 0
        
        lr  = 0.0e0_DP
        i_rho = 1
                  
        do      
            if(i_rho<100) then
                lr  = lr+1.0e-4_DP 
            else if(i_rho<200) then  
                lr  = lr+1.0e-3_DP 
            else if(i_rho<300) then  
                lr  = lr+1.0e-2_DP   
            else if(i_rho<400) then  
                lr  = lr+1.0e-1_DP 
            else if(i_rho<500) then  
                lr  = lr+1.0e0_DP  
            else if(i_rho<1000) then  
                lr  = lr+5.0e0_DP 
            else if(i_rho<1500) then  
                lr  = lr+10.0e0_DP 
            else if(i_rho<2000) then  
                lr  = lr+50.0e0_DP 
            else if(i_rho<2500) then  
                lr  = lr+100.0e0_DP 
            else                             
                go to 30
            end if
                             
            lv2 = 1.0e0_DP/lr
            if(lv2<=1.0e-5_DP) go to 30
            v   = lv2   
            lp2 = Press()    

            if(i_rho>0) then

                !sensible p and within desired p
                if(((lp1<p).and.(lp2>p)).or.((lp1>p).and.(lp2<p))) then
                if((dabs(lp1)<1.0e200_DP).and.(dabs(lp2)<1.0e200_DP))then
                !check g is real
                g_local = (A_ideal( ) + A_mono( ) + A_chain( ) +  A_assoc( ) + A_ion( )) * NA * KB * t + lv2 * lp2 
                if(g_local==g_local) then
                    !interpolate                   
                    grad_local = (lp2 - lp1) / (lv2 - lv1)  
                    v = (p - lp1) / grad_local + lv1
                
                    !Converge the volume   
                    piter=0
                                                       
                    do              
                        v_temp2 = v                                  
                        p_temp2 = Press( )
    
                        if(piter<=plimit1) then                  
                            if(dabs((p_temp2 - p) / p)<P_crit) go to 20
                        else
                            if(dabs((p_temp2 - p) / p)<P_crit2) go to 20             
                        end if
                        
                        d_v = 1.0e-6_DP * v
                        
                        v_temp1 = v - d_v 
                        v = v_temp1
                        p_temp1 = Press( )
                          
                        !interpolate v
                        v = (p-p_temp1) * (v_temp2-v_temp1) / (p_temp2-p_temp1) + v_temp1
    
                        !escape if spurious volume
                        if(v<0.0e0_DP) then 
                            go to 10
                        else if(v/=v) then
                            go to 10
                        end if
                        
                        !escape if too many iterations
                        if(piter>plimit2) go to 20
                        piter = piter + 1                   
                    end do
                    
 20                 nroots = nroots+1
                    lv_array(nroots) = v
                    !lg_array(nroots) = (A_ideal( ) + A_mono( ) + A_chain( ) +  A_assoc( ) + A_ion( )) * NA * KB * t + lv2 * lp2     
                    !better g calc:
                    lg_array(nroots) = 0      
                    do imusum=1,nctypes
                        lg_array(nroots) = lg_array(nroots) + Comp_array(imusum)%xi * Mu(imusum)
                    end do               
                    
                end if
                end if
                end if
            end if
         
 10         lv1 = lv2
            lp1 = lp2
            
            i_rho = i_rho + 1
            
            if(present(limited)) then
                if (i_rho>i_rho_limit) then
                    limited=.false.
                    go to 30
                end if
            end if
            
            if((i_rho>1000).and.(lp2>5e8_DP)) then
                go to 30
            end if
        end do

        !determine most stable volume
 30     if(nroots<1) stop "no volume roots?"
        froot=1        

        if(nroots>1) then
            do i_rho = 2, nroots
                if(lg_array(i_rho)<lg_array(froot)) froot = i_rho
            end do
        end if
       
        v_out = lv_array(froot)

 99     return
    end function Vol_dens_g
!***************************************************************************************************
    function Vol_dens_g2(v_in1, v_in2) result(v_out)  
        implicit none
        
        real(kind=DP), intent(in)   ::  v_in1, v_in2
        real(kind=DP)       ::  r_low, r_high
        real(kind=DP)       ::  v_out
        real(kind=DP)       ::  lv_array(1:10), lg_array(1:10)
        real(kind=DP)       ::  v_temp1, v_temp2, p_temp1, p_temp2, d_v
        real(kind=DP)       ::  lv1, lv2, lp1, lp2, lr, g_local, grad_local
        integer             ::  i_rho, nroots, piter, froot, imusum
        integer, parameter          :: plimit1=100, plimit2=500
        real(kind=DP), parameter    ::  P_crit=1.0e-6_DP,P_crit2=1.0e-4

        nroots = 0
        
        if(1.0e0_DP/v_in1>1.0e0_DP/v_in2) then
            r_high = 2.0e0_DP/v_in1
            r_low  = 0.5e0_DP/v_in2   
        else 
            r_high = 2.0e0_DP/v_in2
            r_low  = 0.5e0_DP/v_in1
        end if
            
        lr  = r_low
        i_rho = 1
                  
        do      
            if(i_rho==1) then
                continue
            else if(i_rho<100) then
                lr  = lr+1.0e-4_DP 
            else if(i_rho<200) then  
                lr  = lr+1.0e-3_DP 
            else if(i_rho<300) then  
                lr  = lr+1.0e-2_DP   
            else if(i_rho<400) then  
                lr  = lr+1.0e-1_DP 
            else if(i_rho<500) then  
                lr  = lr+1.0e0_DP  
            else if(i_rho<1000) then  
                lr  = lr+5.0e0_DP 
            else if(i_rho<1500) then  
                lr  = lr+10.0e0_DP 
            else if(i_rho<2000) then  
                lr  = lr+50.0e0_DP 
            else if(i_rho<2500) then  
                lr  = lr+100.0e0_DP 
            else                             
                go to 30
            end if
            
            if(lr>r_high) go to 30
                            
            lv2 = 1.0e0_DP/lr
            if(lv2<=1.0e-5_DP) go to 30
            v   = lv2   
            lp2 = Press()    

            if(i_rho>0) then

                !sensible p and within desired p
                if(((lp1<p).and.(lp2>p)).or.((lp1>p).and.(lp2<p))) then
                if((dabs(lp1)<1.0e200_DP).and.(dabs(lp2)<1.0e200_DP))then
                !check g is real
                g_local = (A_ideal( ) + A_mono( ) + A_chain( ) +  A_assoc( ) + A_ion( )) * NA * KB * t + lv2 * lp2 
                if(g_local==g_local) then
                    !interpolate                   
                    grad_local = (lp2 - lp1) / (lv2 - lv1)  
                    v = (p - lp1) / grad_local + lv1
                
                    !Converge the volume   
                    piter=0
                                                       
                    do              
                        v_temp2 = v                                  
                        p_temp2 = Press( )
    
                        if(piter<=plimit1) then                  
                            if(dabs((p_temp2 - p) / p)<P_crit) go to 20
                        else
                            if(dabs((p_temp2 - p) / p)<P_crit2) go to 20             
                        end if
                        
                        d_v = 1.0e-6_DP * v
                        
                        v_temp1 = v - d_v 
                        v = v_temp1
                        p_temp1 = Press( )
                          
                        !interpolate v
                        v = (p-p_temp1) * (v_temp2-v_temp1) / (p_temp2-p_temp1) + v_temp1
    
                        !escape if spurious volume
                        if(v<0.0e0_DP) then 
                            go to 10
                        else if(v/=v) then
                            go to 10
                        end if
                        
                        !escape if too many iterations
                        if(piter>plimit2) go to 20
                        piter = piter + 1                   
                    end do
                    
 20                 nroots = nroots+1
                    lv_array(nroots) = v
                    !lg_array(nroots) = (A_ideal( ) + A_mono( ) + A_chain( ) +  A_assoc( ) + A_ion( )) * NA * KB * t + lv2 * lp2     
                    !better g calc:
                    lg_array(nroots) = 0      
                    do imusum=1,nctypes
                        lg_array(nroots) = lg_array(nroots) + Comp_array(imusum)%xi * Mu(imusum)
                    end do               
                    
                end if
                end if
                end if
            end if
         
 10         lv1 = lv2
            lp1 = lp2
            
            i_rho = i_rho + 1
        end do

        !determine most stable volume
 30     if(nroots<1) then
            v_out=1.e10
        else
            froot=1
            
            if(nroots>1) then
                do i_rho = 2, nroots
                    if(lg_array(i_rho)<lg_array(froot)) froot = i_rho
                end do
            end if
            
            v_out = lv_array(froot)
        end if

        return
    end function Vol_dens_g2
!***************************************************************************************************
!***************************************************************************************************
end module Vol_mod
!***************************************************************************************************
!***************************************************************************************************
!
!   Module to add a quote - please add more as and when used
!
!***************************************************************************************************
module Quote_mod
!***************************************************************************************************
contains
!***************************************************************************************************
    subroutine Quote( )
!***************************************************************************************************
!Modules
!=======
        use Zig 
!***************************************************************************************************
    !Variables
    !=========
        implicit none
        
        integer     ::  n_quote, i_quote
        integer     ::  time(1:8), seed
!***************************************************************************************************    
        n_quote = 108  
        time(:)=0        
        
        call DATE_AND_TIME(values=time)     ! Get the current time
        seed = time(1) * (360000*time(5) + 6000*time(6) + 100*time(7) + time(8)) 
  
        call zigset(seed)
        i_quote = anint((n_quote-1)*uni())+1 

        if(i_quote==1) then
            print*, "I am Guybrush Threepwood, mighty pirate. (Monkey Island)"
        else if(i_quote==2) then    
            print*, "You fight like a dairy farmer. (Monkey Island)"
        else if(i_quote==3) then    
            print*, "Look behind you, a three-headed monkey!. (Monkey Island)"
        else if(i_quote==4) then    
            print*, "In my experience, there is no such thing as luck.  "
            print*, " (Obi-Wan Kenobi, Star Wars Episode IV: A New Hope)"
        else if(i_quote==5) then    
            print*, "The needs of the many outweigh the needs of the few. "
            print*, " (Spock, Star Trek II: The Wrath of Khan)"
        else if(i_quote==6) then    
            print*, "You are what you do. A man is defined by his actions, not his memory. "
            print*, " (Kuato, Total Recall)"
        else if(i_quote==7) then    
            print*, "Do what I do. Hold tight and pretend it's a plan! (The Doctor, Doctor Who)"
        else if(i_quote==8) then    
            print*, "Don't panic! (The Hitchhiker's Guide to the Galaxy)"
        else if(i_quote==9) then    
            print*, "Frak! (Battlestar Galactica)"
        else if(i_quote==10) then    
            print*, "Yeah, I got a plan B: making plan A work! (Stingray)"
        else if(i_quote==11) then    
            print*, "Why don't you smegging well smeg off you annoying little smeggy smegging smegger! "
            print*, " (Red Dwarf)"
        else if(i_quote==12) then    
            print*, "Make it so. (Star Trek TNG)"
        else if(i_quote==13) then    
            print*, "Hokey religions and ancient weapons are no match for a good blaster at your side, kid. "
            print*, " (Hans Solo, Star Wars Episode IV: A New Hope)"
        else if(i_quote==14) then    
            print*, "Roads? Where we're going, we don't need roads. (Back to the Future)"
        else if(i_quote==15) then    
            print*, "Come with me if you want to live. (The Terminator)"
        else if(i_quote==16) then    
            print*, "We'd better get back, 'cause it'll be dark soon, and they mostly come at night... mostly. "
            print*, " (Aliens)"
        else if(i_quote==17) then    
            print*, "Dead or alive, you're coming with me! (Robocop)"
        else if(i_quote==18) then    
            print*, "I mean, have you ever actually seen your brain? (Ghost in the Shell)"
        else if(i_quote==19) then    
            print*, "We spend $250 billion a year on defence, and here we are! The fate of the planet "
            print*, " is in the hands of a bunch of retards I wouldn't trust with a potato gun! (Armageddon)"
        else if(i_quote==20) then    
            print*, "Those who know do not speak. Those who speak do not know. (Lao Tzu, Tao Teh Ching)" 
        else if(i_quote==21) then    
            print*, "The truth is not always beautiful, nor beautiful words the truth. (Lao Tzu, Tao Teh Ching)" 
        else if(i_quote==22) then    
            print*, "When I let go of what I am, I become what I might be. (Lao Tzu, Tao Teh Ching)" 
        else if(i_quote==23) then    
            print*, "Care about what other people think and you will always be their prisoner. "
            print*, " (Lao Tzu, Tao Teh Ching)" 
        else if(i_quote==24) then    
            print*, "A man with outward courage dares to die; a man with inner courage dares to live. "
            print*, " (Lao Tzu, Tao Teh Ching)" 
        else if(i_quote==25) then    
            print*, "Teenage angst has paid off well, Now I'm bored and old. (Nirvana)"
        else if(i_quote==26) then    
            print*, "Hey! Wait! I got a new complaint. (Nirvana)"
        else if(i_quote==27) then    
            print*, "I don't care what you think unless it's about me. (Nirvana)"
        else if(i_quote==28) then    
            print*, "Just because you're paranoid. Don't mean they're now after you. (Nirvana)"
        else if(i_quote==29) then    
            print*, "I'm so happy cause today I've found my friends. They're in my head. (Nirvana)"
        else if(i_quote==30) then    
            print*, "Please allow me to introduce myself / I'm a man of wealth and taste/ I've been"
            print*, " around for a long, long year / Stole many a man's soul and faith. (The Rolling Stones)"
        else if(i_quote==31) then    
            print*, "There's no time for us. There's no place for us. What is this thing that builds our "
            print*, " dreams yet slips away from us? (Queen)"
        else if(i_quote==32) then    
            print*, "And when the brokenhearted people living in the world agree, there will be an answer,"
            print*, " let it be. (The Beatles)"
        else if(i_quote==33) then    
            print*, "The fewer the facts, the stronger the opinion. (Arnold H. Glasow)"
        else if(i_quote==34) then    
            print*, "The distance between insanity and genius is measured only by success. (Bruce Feirstein)"
        else if(i_quote==35) then    
            print*, "Your theory is crazy, but it's not crazy enough to be true. (Niels Bohr)"
        else if(i_quote==36) then    
            print*, "Science is nothing but perception. (Plato)"
        else if(i_quote==37) then    
            print*, "If your experiment needs statistics, you ought to have done a better experiment. "
            print*, " (Ernest Rutherford)"
        else if(i_quote==38) then    
            print*, "A physicist is an atom's way of knowing about atoms. (George Wald)"
        else if(i_quote==39) then    
            print*, "Valid criticism does you a favour. (Carl Sagan)"
        else if(i_quote==40) then    
            print*, "Chemistry, unlike other sciences, sprang originally from delusions and superstitions, "
            print*, " and was at its commencement exactly on a par with magic and astrology. (Thomas Thomson)"
        else if(i_quote==41) then    
            print*, "To the electron -- may it never be of any use to anybody. (JJ. Thomson)"
        else if(i_quote==42) then    
            print*, "After climbing a great hill, one only finds there are many more hills to climb. "
            print*, " (Nelson Mandela)"
        else if(i_quote==43) then    
            print*, "If you cannot do great things, do small things in a great way. (Napoleon Hill)"
        else if(i_quote==44) then    
            print*, "It is a rough road that leads to the heights of greatness. (Lucius Annaeus Seneca)"
        else if(i_quote==45) then    
            print*, "Only those who dare to fail greatly can ever achieve greatly. (Robert Kennedy)"
        else if(i_quote==46) then    
            print*, "It is hard to be humble when you're as great as I am. (Muhammad Ali)"
        else if(i_quote==47) then    
            print*, "Two roads diverged in a wood, and I - I took the one less traveled by, and that "
            print*, " has made all the difference. (Robert Frost)"
        else if(i_quote==48) then    
            print*, "Life is what happens to you while you're busy making other plans. (John Lennon)"
        else if(i_quote==49) then    
            print*, "The most common way people give up their power is by thinking they don't have any."
            print*, " (Alice Walker)"
        else if(i_quote==50) then    
            print*, "The best time to plant a tree was 20 years ago. The second best time is now. "
            print*, " (Chinese Proverb)"
        else if(i_quote==51) then    
            print*, "Every child is an artist.  The problem is how to remain an artist once he grows up. "
            print*, " (Pablo Picasso)"
        else if(i_quote==52) then    
            print*, "Whether you think you can or you think you can't, you're right. (Henry Ford)"
        else if(i_quote==53) then    
            print*, "People often say that motivation doesn't last. Well, neither does bathing.  "
            print*, " That's why we recommend it daily. (Zig Ziglar)"
        else if(i_quote==54) then    
            print*, "There is only one way to avoid criticism: do nothing, say nothing, and be nothing. "
            print*, " (Aristotle)"
        else if(i_quote==55) then    
            print*, "Everything you've ever wanted is on the other side of fear. (George Addair)"
        else if(i_quote==56) then    
            print*, "Teach thy tongue to say, I do not know, and thou shalt progress. (Maimonides)"
        else if(i_quote==57) then    
            print*, "How wonderful it is that nobody need wait a single moment before starting to "
            print*, " improve the world. (Anne Frank)"
        else if(i_quote==58) then    
            print*, "If the wind will not serve, take to the oars. (Latin Proverb)"
        else if(i_quote==59) then    
            print*, "Challenges are what make life interesting and overcoming them is what makes life"
            print*, " meaningful. (Joshua J. Marine)"
        else if(i_quote==60) then    
            print*, "The person who says it cannot be done should not interrupt the person who is "
            print*, " doing it. (Chinese Proverb)"
        else if(i_quote==61) then    
            print*, "Build your own dreams, or someone else will hire you to build theirs. (Farrah Gray)"
        else if(i_quote==62) then    
            print*, "It does not matter how slowly you go as long as you do not stop. (Confucius)"
        else if(i_quote==63) then    
            print*, "Everything in moderation, including moderation. (Oscar Wilde)"
        else if(i_quote==64) then    
            print*, "Just taught my kids about taxes by eating 38% of their ice cream. (Conan O'Brien)"
        else if(i_quote==65) then    
            print*, "I dream of a better tomorrow, where chickens can cross the road and not be questioned"
            print*, " about their motives. (Anonymous)"
        else if(i_quote==66) then    
            print*, "It is the mark of an educated mind to be able to entertain a thought without "
            print*, " accepting it. (Aristotle)"
        else if(i_quote==67) then    
            print*, "You see things and you say Why? But I dream things that never were and I say "
            print*, " Why not? (George Bernard Shaw)"
        else if(i_quote==68) then    
            print*, "I have never in my life learned anything from a man who agreed with me. "
            print*, " (Dudley Field Malone)"
        else if(i_quote==69) then    
            print*, "We have all heard that a million monkeys banging on a million typewriters will "
            print*, " eventually reproduce the entire works of Shakespeare. Now, thanks to the internet, "
            print*, " we know this is not true. (Robert Wilensky)"
        else if(i_quote==70) then    
            print*, "Every truth passes through three stages before it is recognised. In the first,"
            print*, " it is ridiculed. In the second, it is opposed. In the third it is regarded as "
            print*, " self-evident. (Arthur Schopenhauer)"
        else if(i_quote==71) then    
            print*, "I always arrive late at the office, but I make up for it by leaving early. "
            print*, " (Charles Lamb)"
        else if(i_quote==72) then    
            print*, "The elevator to success is out of order. You'll have to use the stairs... "
            print*, " one step at a time. (Joe Girard)"
        else if(i_quote==73) then    
            print*, "The greatest way to live with honour in this world is to be what we"
            print*, " pretend to be. (Socrates)"
        else if(i_quote==74) then    
            print*, "Laws control the lesser man... Right conduct controls the greater one. (Mark Twain)"
        else if(i_quote==75) then    
            print*, "Great men are seldom over-scrupulous in the arrangement of their attire. "
            print*, " (Charles Dickens)"
        else if(i_quote==76) then    
            print*, "Time's fun when you're having flies. (Kermit the Frog )" 
        else if(i_quote==77) then    
            print*, "You have to be odd to be number one. (Dr. Seuss)"
        else if(i_quote==78) then    
            print*, "You have brains in your head. You have feet in your shoes. You can steer "
            print*, " yourself any direction you choose. (Dr. Seuss)"
        else if(i_quote==79) then    
            print*, "Sometimes the questions are complicated and the answers are simple. (Dr. Seuss)"
        else if(i_quote==80) then    
            print*, "Today you are you. That is truer than true. There is no one alive who is "
            print*, " youer than you. (Dr. Seuss)"
        else if(i_quote==81) then    
            print*, "A little nonsense now and then, is relished by the wisest men. (Roald Dahl)"
        else if(i_quote==82) then    
            print*, "In ancient times cats were worshipped as gods. They have not forgotten this. "
            print*, "(Terry Pratchett)"
        else if(i_quote==83) then    
            print*, "The trouble with having an open mind, of course, is that people will insist on"
            print*, " coming along and trying to put things in it. (Terry Pratchett)"
        else if(i_quote==84) then    
            print*, "The pen is mightier than the sword if the sword is very short, and the pen is "
            print*, " very sharp. (Terry Pratchett)"
        else if(i_quote==85) then    
            print*, "The presence of those seeking the truth is infinitely to be preferred to the "
            print*, " presence of those who think they've found it. (Terry Pratchett)"
        else if(i_quote==86) then      
            print*, " We have created amazing machines... and used them to destroy our fellow men. "
            print*, " (Metropolis)"
        else if(i_quote==87) then      
            print*, "You're mad. Bonkers. Off your head... but I'll tell you a secret... "
            print*, "all of the best people are. (Alice in Wonderland)"
        else if(i_quote==88) then      
            print*, "Don't cry because it's over, smile because it happened. (Dr. Seuss)"
        else if(i_quote==89) then      
            print*, "We are all in the gutter, but some of us are looking at the stars."
            print*, "(Oscar Wilde)"
        else if(i_quote==90) then      
            print*, "I may not have gone where I intended to go, but I think I have ended"
            print*, "up where I needed to be. (Douglas Adams)"
        else if(i_quote==91) then      
            print*, "For every minute you are angry you lose sixty seconds of happiness."
            print*, "(Ralph Walso Emerson)"
        else if(i_quote==92) then      
            print*, "If you can't explain it to a six year old, you don't understand it "
            print*, "yourself. (Albert Einstein)"
        else if(i_quote==93) then      
            print*, "Success is not final, failure is not fatal: it is the courage to "
            print*, "continue that counts. (Winston Churchill)"
        else if(i_quote==94) then      
            print*, "Reality continues to ruin my life. (Bill Watterson)"
        else if(i_quote==95) then      
            print*, "I am free of all prejudice. I hate everyone equally. (W.C. Fields)"
        else if(i_quote==96) then      
            print*, "It's no use going back to yesterday, because I was a different"
            print*, "person then. (Lewis Carroll)"
        else if(i_quote==97) then      
            print*, "What you're supposed to do when you don't like a thing is change"
            print*, "it. If you can't change it, change the way you think about it. "
            print*, "Don't complain. (Maya Angelou)"
        else if(i_quote==98) then 
            print*, "Experience is that marvelous thing that enables you to recognize "
            print*, "a mistake when you make it again. (F. P. Jones)"
        else if(i_quote==99) then 
            print*, "Human beings, who are almost unique in having the ability to learn "
            print*, "from the experience of others, are also remarkable for their apparent "
            print*, "disinclination to do so. (Douglas Adams)"
        else if(i_quote==100) then
            print*, "My opinions may have changed, but not the fact that I am right. (Ashleigh Brilliant)"
        else if(i_quote==101) then
            print*, "Every man is born as many men and dies as a single one. (Martin Heidegger)"
        else if(i_quote==102) then
            print*, "Many are stubborn in pursuit of the path they have chosen,"
            print*, "few in pursuit of the goal. (Friedrich Nietzsche)"
        else if(i_quote==103) then
            print*, "If you are lonely when you're alone, you are in bad company. (Jean-Paul Sartre)"
        else if(i_quote==104) then
            print*, "There are two things a person should never be angry at, what they can help,"
            print*, "and what they cannot. (Plato)"
        else if(i_quote==105) then
            print*, "No tree, it is said, can grow to heaven unless its roots reach down to hell. (Carl Jung)"
        else if(i_quote==106) then
            print*, "Truth often suffers more by the heat of its defenders than "
            print*, "the arguments of its opposers. (William Penn)"
        else if(i_quote==107) then
            print*, "You can't depend on your eyes when your imagination is out of focus. (Mark Twain)"
        else if(i_quote==108) then
            print*, "Sometimes you face difficulties not because you're doing something wrong,"
            print*, "but because you're doing something right. (Joel Osteen)"
        end if    
        
        return
    end subroutine Quote
end module Quote_mod

!***************************************************************************************************

    Module c05qbfe_mod

!     C05QBF Example Program Module:
!            Parameters and User-defined Routines

!     .. Use Statements ..
      Use nag_library, Only: nag_wp
	  use Types           ! Definitions of types and double precision
	  use Global          ! Important global parameters
	  use Pressure
	  use ChemPot
    use Vol_mod
        !     .. Implicit None Statement ..
      Implicit None
!     .. Accessibility Statements ..
      Private
      Public                           :: fcn
!     .. Parameters ..
      Integer, Parameter, Public       :: n_vle = 3, nout_vle = 6
	  Real (Kind=nag_wp), Public			   :: P_L, P_V, Mu_L_1, Mu_L_2, Mu_V_1, Mu_V_2, x_1, x_init_1, x_init_2, x_init_3, x_init_4
    Contains
      Subroutine fcn(n_vle,x,fvec,iuser,ruser,iflag)

!       .. Scalar Arguments ..
        Integer, Intent (Inout)        :: iflag
        Integer, Intent (In)           :: n_vle
!       .. Array Arguments ..
        Real (Kind=nag_wp), Intent (Out) :: fvec(n_vle)
        Real (Kind=nag_wp), Intent (Inout) :: ruser(*)
        Real (Kind=nag_wp), Intent (In) :: x(n_vle)
        Integer, Intent (Inout)        :: iuser(*)
!       .. Executable Statements ..
		v   = x(1)				! Liquid volume is x1
		Comp_array(1)%xi = x_1 !x1
		Comp_array(2)%xi = 1.0_nag_wp - x_1 !x2
		p = Press()
		P_L = p
		Mu_L_1 = Mu(1)
		Mu_L_2 = Mu(2)
		v   = x(2)				! Vapor volume is x2
		Comp_array(1)%xi = x(3) ! y1
		Comp_array(2)%xi =1.0_nag_wp - x(3) !y2
		p = Press()
		P_V = p
		Mu_V_1 = Mu(1)
		Mu_V_2 = Mu(2)
		
		
		
        fvec(1) = P_L - P_V
        fvec(2) = Mu_L_1 - Mu_V_1
		fvec(3) = Mu_L_2 - Mu_V_2
		
        !       Set iflag negative to terminate execution for any reason.
        iflag = 0
        Return
      End Subroutine fcn
    End Module c05qbfe_mod
!***************************************************************************************************    
subroutine c05qbfe(output_vle)

!     C05QBF Example Main Program

!     .. Use Statements ..
      Use c05qbfe_mod, Only: fcn, n_vle, nout_vle, x_init_1, x_init_2, x_init_3, x_init_4

      Use nag_library, Only: c05qbf, dnrm2, nag_wp, x02ajf
!     .. Implicit None Statement ..
      Implicit None
!     .. Local Scalars ..
      Real (Kind=nag_wp)               :: fnorm, xtol
      Integer                          :: i, ifail
!     .. Local Arrays ..
      Real (Kind=nag_wp), Allocatable  :: fvec(:), x(:)
      Real (Kind=nag_wp)               :: ruser(1)
      Integer                          :: iuser(1)
      !     .. Intrinsic Procedures ..
      Intrinsic                        :: sqrt
	  Real  (Kind=nag_wp)              :: output_vle(3)
!     .. Executable Statements ..

      Allocate (fvec(n_vle),x(n_vle))

!     The following starting values provide a rough solution.

      x(1) = x_init_1
      x(2) = x_init_2
      x(3) = x_init_3
      
      xtol = sqrt(x02ajf())

      ifail = -1
      Call c05qbf(fcn,n_vle,x,fvec,xtol,iuser,ruser,ifail)

      If (ifail==0 .Or. ifail==2 .Or. ifail==3 .Or. ifail==4) Then
        If (ifail==0) Then
!         The NAG name equivalent of dnrm2 is f06ejf
          fnorm = dnrm2(n_vle,fvec,1)
        Else
          Write (nout_vle,*)
          Write (nout_vle,*) 'Approximate solution'
		  x(3)=10
        End If
        !Write (nout,99998)(x(i),i=1,n)
		Do i = 1, n_vle
           output_vle(i) = x(i)
        End Do
      End If

99999 Format (1X,A,E12.4)
99998 Format (1X,3E12.4)
  return
End

!   E04JCF Example Program Text
!   Mark 26.1 Release. NAG Copyright 2017.
    Module e04jcfe_mod

!     E04JCF Example Program Module:
!            Parameters and User-defined Routines

!     .. Use Statements ..
      Use nag_library, Only: nag_wp
	  use Types           ! Definitions of types and double precision
	  use Global          ! Important global parameters
	  use Pressure
	  use Input
	  use Input_opt   ! Read optimisation parameters
	  use Vol_mod
	  use c05qbfe_mod
	  use ChemPot
!     .. Implicit None Statement ..
      Implicit None
!     .. Accessibility Statements ..
      Private
      Public                           :: monfun, objfun
!     .. Parameters ..
      Integer, Parameter, Public       :: nout = 6 
      ! integer, public                         ::  n_opt, opt_max
      !real(kind=DP),allocatable, public       ::  min_bound(:), max_bound(:)
      
    Contains
      Subroutine objfun(n,x,f,iuser,ruser,inform)

!       .. Parameters ..
        ! Real (Kind=nag_wp), Parameter  :: five = 5.0_nag_wp
        ! Real (Kind=nag_wp), Parameter  :: ten = 1.0E1_nag_wp
        ! Real (Kind=nag_wp), Parameter  :: two = 2.0_nag_wp
!       .. Scalar Arguments ..
        Real (Kind=nag_wp), Intent (Out) :: f
        Integer, Intent (Out)          :: inform
        Integer, Intent (In)           :: n
!       .. Array Arguments ..
        Real (Kind=nag_wp), Intent (Inout) :: ruser(*)
        Real (Kind=nag_wp), Intent (In) :: x(n)
        Integer, Intent (Inout)        :: iuser(*)
		integer 						:: i
		real(Kind=nag_wp)				::  p_oj, y_oj(2),values(3), mu_1, mu_2, &
										& p_squ, y_squ(2)
!       .. Executable Statements ..
        inform = 0
		
		!allocate(min_bound(1:opt_num),max_bound(1:opt_num))
		do i = 1, opt_num
			if(param_key(i)==1) then               
				!min_bound(i)=min_num(i) * ANG
				!max_bound(i)=max_num(i) * ANG
				sig(param_index(i,1),param_index(i,1)) = x(i)* ANG/1000  !1000 is to make nag get converged
        !print*, Seg_array(param_index(i,1))%sig, x(i)
            else if(param_key(i)==2) then
				!min_bound(i)=min_num(i)
				!max_bound(i)=max_num(i)
				eps(param_index(i,1),param_index(i,1)) = x(i)/10     !10 is to make nag get converged 
            else if(param_key(i)==3) then
				!min_bound(i)=min_num(i)
				!max_bound(i)=max_num(i)
				lr(param_index(i,1),param_index(i,1)) = x(i)/100   !100 is to make nag get converged
            else if(param_key(i)==4) then
				!min_bound(i)=min_num(i)
				!max_bound(i)=max_num(i)
				la(param_index(i,1),param_index(i,1)) = x(i)
            else if(param_key(i)==5) then
				!min_bound(i)=min_num(i)
				!max_bound(i)=max_num(i)
				Seg_array(param_index(i,1))%sf = x(i)/10000 !10000 is to make nag get converged
            else if(param_key(i)==6) then
				!min_bound(i)=min_num(i)
				!max_bound(i)=max_num(i)
				Seg_array(param_index(i,1))%nseg = x(i)
			else if(param_key(i)==7)then
				!min_bound(i)=min_num(i)
				!max_bound(i)=max_num(i)
				eps(param_index(i,1), param_index(i,2))=x(i)/10  !10 is to make nag get converged
				eps(param_index(i,2), param_index(i,1))=x(i)/10 !10 is to make nag get converged				
            else if(param_key(i)==8)then
				!min_bound(i)=min_num(i)
				!max_bound(i)=max_num(i)
				lr(param_index(i,1), param_index(i,2))=x(i)/100 !100 is to make nag get converged
				lr(param_index(i,2), param_index(i,1))=x(i)/100 !100 is to make nag get converged
            else if(param_key(i)==9)then
				!min_bound(i)=min_num(i)
				!max_bound(i)=max_num(i)
				ehb(param_index(i,1),param_index(i,2), &
				& param_index2(i,1),param_index2(i,2)) = x(i)
				ehb(param_index2(i,1),param_index2(i,2), &
				& param_index(i,1),param_index(i,2)) = x(i)
            else if(param_key(i)==10)then
				!min_bound(i)=min_num(i) * ANG**3.0_DP * NA
				!max_bound(i)=max_num(i) * ANG**3.0_DP * NA
				khb(param_index(i,1),param_index(i,2), &
				& param_index2(i,1),param_index2(i,2)) = x(i)* ANG**3.0_DP * NA/10 !10 is to make nag get converged
				khb(param_index2(i,1),param_index2(i,2), &
				& param_index(i,1),param_index(i,2)) = x(i)* ANG**3.0_DP * NA/10 !10 is to make nag get converged
			end if
		end do
		
		p_oj = 0
		y_oj(1)=0
		y_oj(2)=0
		p_squ=0
		y_squ(1)=0
		y_squ(2)=0
		


	if (properties%opt_l(6)) then
		do i = 1, properties%nmu 
          Comp_array(1:nctypes)%xi = properties%ximu(i,1:nctypes)
		  t = properties%t_opt(i)
		  p = properties%p_opt(i)
          x_1 = properties%ximu(i,1)
		  P = 2*P
		  v = Vol_dens_g( )
		  mu_1 = Mu_res(1)
		  mu_2 = Mu_res(2)
		  x_init_1 = v
		  x_init_2 = 8.314*t/(properties%p_opt(i))
		  x_init_3 = exp(mu_1/(8.314*t))*x_1/(exp(mu_1/(8.314*t))*x_1+exp(mu_2/(8.314*t))*(1-x_1))
		  p = properties%p_opt(i)

          call c05qbfe(values )
			p_oj=p_oj+ABS(((P_L+P_V)/2-properties%p_opt(i))/properties%p_opt(i))
			p_squ=p_squ+(((P_L+P_V)/2-properties%p_opt(i))/properties%p_opt(i))**2
			y_oj(1)=y_oj(1)+ABS((values(3)-properties%yimu(i,1))/properties%yimu(i,1))
			y_squ(1)=y_squ(1)+((values(3)-properties%yimu(i,1))/properties%yimu(i,1))**2
			y_oj(2)=y_oj(2)+ABS((1.0_nag_wp-values(3)-properties%yimu(i,2))/properties%yimu(i,2))
			y_squ(2)=y_squ(2)+((1.0_nag_wp-values(3)-properties%yimu(i,2))/properties%yimu(i,2))**2
		end do
			p_oj=p_oj/properties%nmu
			p_squ=p_squ/properties%nmu
			y_oj(1)=y_oj(1)/properties%nmu
			y_squ(1)=y_squ(1)/properties%nmu
			y_oj(2)=y_oj(2)/properties%nmu
			y_squ(2)=y_squ(2)/properties%nmu
		end if		
        
		
		f = p_squ + y_squ(1) + y_squ(2)
   !print*, "f          x(i)"
   write(*, 101)"f=",f, "p_oj=", p_oj, "y_oj1=", y_oj(1), "y_oj2=", y_oj(2),"x(i)=",x(1:opt_num)
101 Format ((2X,A2),(2X,F6.4),(2X,A5),(2X,F6.4),(2X,A6),(2X,F6.4),(2X,A6),(2X,F6.4),(2X,A5),*(2X,F10.4))

        Return

      End Subroutine objfun
      Subroutine monfun(n,nf,x,f,rho,iuser,ruser,inform)

!       .. Scalar Arguments ..
		!use Types       ! Definitions of types and double precision
		!use Global,only: opt_num, min_num, max_num, init_values      ! Important global parameters
		
        Real (Kind=nag_wp), Intent (In) :: f, rho
        Integer, Intent (Out)          :: inform
        Integer, Intent (In)           :: n, nf
!       .. Array Arguments ..
        Real (Kind=nag_wp), Intent (Inout) :: ruser(*)
        Real (Kind=nag_wp), Intent (In) :: x(n)
        Integer, Intent (Inout)        :: iuser(*)
!       .. Local Scalars ..
        Logical                        :: verbose_output
!       .. Executable Statements ..
        inform = 0

        Write (nout,Fmt=99999) 'Monitoring: new trust region radius =', rho

!       Set this to .True. to get more detailed output
        verbose_output = .false.

        If (verbose_output) Then
          Write (nout,Fmt=99998) 'Number of function calls =', nf
          Write (nout,Fmt=99997) 'Current function value =', f
          Write (nout,Fmt=99996) 'The corresponding X is:', x(1:n)
        End If

        Return
99999   Format (/,4X,A,1P,E13.3)
99998   Format (4X,A,I16)
99997   Format (4X,A,1P,E12.4)
99996   Format (4X,A,/,(4X,5E12.4))
      End Subroutine monfun
    End Module e04jcfe_mod

module nag_e04jcfe_mod
contains    
	subroutine e04jcfe(output)

!     Example problem for E04JCF.

!     .. Use Statements ..
      Use e04jcfe_mod, Only: monfun, nout, objfun!, min_bound, max_bound
      Use nag_library, Only: e04jcf, nag_wp, x02alf
	  !use Types       ! Definitions of types and double precision
	  use Global,only: opt_num, min_num, max_num, init_values, param_key, ANG, NA       ! Important global parameters
!     .. Implicit None Statement ..
      Implicit None
!     .. Local Scalars ..
      Real (Kind=nag_wp)               :: f, infbnd, rhobeg, rhoend
      Integer                          :: ifail, maxcal, n, nf, npt
!     .. Local Arrays ..
      Real (Kind=nag_wp), Allocatable  :: bl(:), bu(:), x(:), output(:)
      Real (Kind=nag_wp)               :: ruser(1)!, output(opt_num)
      Integer                          :: iuser(1), i
!     .. Executable Statements ..
      Write (nout,*) 'E04JCF Example Program Results'

      maxcal = 5000
      rhobeg = 100
      rhoend = 1.0E-6_nag_wp
      n = opt_num !4
      npt = ((n+1)*(n+2))/2   !2*n + 1

!     x(3) is unconstrained, so we're going to set bl(3) to a large
!     negative number and bu(3) to a large positive number.

      infbnd = x02alf()**0.25_nag_wp
	!print *, infbnd
      Allocate (bl(n),bu(n),x(n))
		do i=1, opt_num
			bl(i) = min_num(i) !(/1.0_nag_wp,-2.0_nag_wp,-infbnd,1.0_nag_wp/)
			bu(i) = max_num(i) !(/3.0_nag_wp,0.0_nag_wp,infbnd,3.0_nag_wp/)
			x(i) = init_values(i) !(/3.0_nag_wp,-1.0_nag_wp,0.0_nag_wp,1.0_nag_wp/)
			!print *, bl(i),bu(i),x(i)
		end do
		!print*, bl(1:n),bu(1:n),x(1:n),output(1:n)
   print*, "initial values"
   write(*,102) x(1:n)
102 Format (*(2X,F9.4))
      ifail = -1
      Call e04jcf(objfun,n,npt,x,bl,bu,rhobeg,rhoend,monfun,maxcal,f,nf,iuser, &
        ruser,ifail)
	!print*, "test"
      Select Case (ifail)
      Case (0,2:5)
	  !print*, "test"
        If (ifail==0) Then
          Write (nout,Fmt=99999) 'Successful exit from E04JCF.',               &
            'Function value at lowest point found =', f
        Else
          Write (nout,Fmt=99998)                                               &
            'On exit from E04JCF, function value at lowest point found =', f
        End If
		!Write (nout,Fmt=99997) 'The corresponding X is:', x(1:n)
		!print*, "test"
		!Allocate (output(n))	
		do i = 1, opt_num
			
			!print*, "test", x(i)
				output(i) = x(i)
			!print*, "test2", x(i)	
		end do
	!print*, "test"
      End Select

99999 Format (2(/,1X,A),1P,E13.3)
99998 Format (/,1X,A,1P,E13.3)
99997 Format (1X,A,/,(2X,5E13.3))
	return
 End
 end 
 
!***************************************************************************************************
!
!       MANSAFT OPTIMISER
!       =================
!
!       Version 2.7
!***************************************************************************************************
!   SAFT gamma mie program to:
!       1. Read general parameter list form input and
!       2. Calculate system parameters
!       3. Optimise parameters according to simplex method
!***************************************************************************************************
!   Created by SJ Halstead Oct 2015
!       Update 2  Jan 2016
!           - tidied and standardised units
!       Update 3  Feb 2016
!           - added simple P / V calculator
!       Update 4  Mar 2016
!           - standardised MONO units
!       Update 5  Mar 2016
!           - reorganised to focus on propety calculation
!       Update 6  May 2016
!           - property calculations enabled for:
!               - P/V (able to select correct V using lowest G)
!               - Chem pot
!               - Phase equilibrium
!       Update 7  Sept 2016
!           - restructured with FORTRAN best practice
!       Update 7
!           - added electrolyte optimiser
!***************************************************************************************************
!       REFERENCES
!       ==========
!
!       [1] Accurate statistical associating fluid theory for chain molecules formed from Mie segments
!               J Chem Phys 193, 154504 (2013)
!               T. Lafitte, A. Apostolakou, C. Avendano, A. Galindo, C. Adjiman, E. Muller, G. Jackson
!
!               (Main reference)
!
!
!       [2] Prediction of thermodynamic propertied and phase behavior of fluids and mixtures with the
!           SAFT-gamma Mie group-contribution equation of state
!               J. Chem. Eng. Data, 59, 3272-3288 (2014)
!               S. Dufal, V. Papaioannou, M. Sadeqzadeh, T> Pogiatzis, A. Chremos, C.S. Adjiman, 
!               G. Jackson, A. Galindo
!
!               (Associations)
!
!       
!       [3] The A in SAFT: developing the contribution of association to the Helmholtz free energy
!           within a Wertheim TPT1 treatment of generic Mie fluids
!               Mol. Phys. 113, 948-984 (2015)
!               S. Dufal, T. Lafitte, A.J. Haslam, A. Galindo, G.N.I. Clark, C.Vega, G. Jackson            
!
!               (Association coefficients)            
!
!
!***************************************************************************************************

program optimiser
!***************************************************************************************************
!Modules
!=======
    use Types       ! Definitions of types and double precision
    use Input       ! Read the input
    use Input_opt   ! Read optimisation parameters
    use Global      ! Important global parameters
    !use Simplex_mod     ! Simplex optimiser
    use Quote_mod       ! Random quote
    use nag_e04jcfe_mod 
	
!***************************************************************************************************
!Variables
!=========
    implicit none
    
    !Clock
    integer             ::  start, finish, rate
    !Property calculation
    character(len=1)    ::  prop
    !Counter
    integer             ::  i
    !Optimiser
    integer                         ::  n_opt !, opt_max
    real(kind=DP),allocatable       ::  min_opt(:), max_opt(:), paramin(:,:), paramout(:), values(:)!, &
									!&	min_init(:), max_init(:)
    real(kind=DP)                   ::  errout!, values(opt_num)
	! real(kind=DP),allocatable,public::  init_values(:)
!***************************************************************************************************
!Program header
!==============
    print*, " "
    print*, "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    print*, "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    print*, "++                                                                             ++"
    print*, "++                                                                             ++"
    print*, "++            M     M    A    N   N   SSSS    A    FFFFF  TTTTT                ++"
    print*, "++            MM   MM   A A   NN  N  SS      A A   F        T                  ++"
    print*, "++            M M M M  A   A  N N N   SSS   A   A  FFF      T                  ++"
    print*, "++            M  M  M  A A A  N  NN     SS  A A A  F        T                  ++"
    print*, "++            M     M  A   A  N   N  SSSS   A   A  F        T                  ++"
    print*, "++                                                                             ++"
    print*, "++                                                                             ++"
    print*, "++                              MANSAFT                                        ++"
    print*, "++                     Manchester SAFT Gamma Mie                               ++"
    print*, "++                                                                             ++"
    print*, "++                                                                             ++"
    print*, "++                      ****     *****     ********                            ++"
    print*, "++                     ******    ******    ********                            ++"
    print*, "++                    **    **   **  ***      **                               ++"
    print*, "++                    **    **   **  ***      **                               ++"
    print*, "++                    **    **   ******       **                               ++"
    print*, "++                    **    **   ***          **                               ++"
    print*, "++                     ******    ***          **                               ++"
    print*, "++                      ****     ***          **                               ++"
    print*, "++                                                                             ++"
    print*, "++                VERSION 2.7 - The Ant on the Elephant                        ++"
    print*, "++                                2018                                         ++"
    print*, "++                                                                             ++"
    print*, "++                                                                             ++"
    print*, "++                       S.J. Halstead & A.J.Masters                           ++"
    print*, "++                                                                             ++"
    print*, "++         School of Chemical Engineering and Analytical Science               ++"
    print*, "++                     The University of Manchester                            ++"
    print*, "++                                                                             ++"
    print*, "++               contact:  simon.halstead@manchester.ac.uk                     ++"
    print*, "++                         simon.halstead@gmail.com                            ++"
    print*, "++                                                                             ++"
    print*, "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    print*, "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    print*, " "
!***************************************************************************************************
!Start clock
!===========
    call System_clock( start, rate )
    print*,"optimiser"    
!***************************************************************************************************
!Read input
!==========
    call Read_input( )    
     
    if((properties%type/='opt').and.(properties%type/='OPT').and.(properties%type/='Opt')   &
    & .and.(properties%type/='eopt').and.(properties%type/='EOPT')                          &
    & .and.(properties%type/='ilopt').and.(properties%type/='ILOPT')) then
        stop "For optimisation, input run type must be OPT, EOPT or ILOPT"
    end if

    call Read_opt( n_opt, max_opt, min_opt, paramin ) 
		!print*, opt_num, min_num(1:opt_num), max_num(1:opt_num)
!***************************************************************************************************    
		!allocate( init_values(1:opt_num))!, min_init(1:opt_num), max_init(1:opt_num))
		
		!do i = 1, opt_num
			!init_values(i) = (min_num(i) + max_num(i))/2
			!print *, init_values(i)
		!end do
	 !print *, init_values(1:opt_num)
	
	allocate(values(1:opt_num))
	call e04jcfe(values)
	do i = 1, opt_num
			if(param_key(i)==1) then
				values(i) = values(i)/1000  !1000 is to make nag get converged
        !print*, Seg_array(param_index(i,1))%sig, x(i)
            else if(param_key(i)==2) then
				values(i) = values(i)/10     !10 is to make nag get converged 
            else if(param_key(i)==3) then
				values(i) = values(i)/100   !100 is to make nag get converged
            else if(param_key(i)==5) then
				values(i) = values(i)/10000 !10000 is to make nag get converged
			else if(param_key(i)==7)then
				values(i) = values(i)/10  !10 is to make nag get converged				
            else if(param_key(i)==8)then
				values(i) = values(i)/100 !100 is to make nag get converged
            else if(param_key(i)==10)then
				values(i) = values(i)/10 !10 is to make nag get converged
			end if
		end do
	write(*,100)  values(1:opt_num)

100 Format (2X,F9.4,$)
!***************************************************************************************************
!End time
!========
    call System_clock(finish)
   
    print*,  " "
    print*,  "##################################################################"
    print*,  "                    ALL DONE  "
    print*,  "    Time taken: ",REAL(finish-start)/real(rate)," seconds"
    print*,  "##################################################################"
    print*,  " "
    print*,  " "
    
    call Quote( )
    print*, " "
    print*, "#############################################################################################"
    print*, " "
    print*, " "
    print*, " "
    print*, " "
    print*, " "

!***************************************************************************************************            
    stop
!***************************************************************************************************

!***************************************************************************************************
end program optimiser
!
