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
        use Zig_mod 
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