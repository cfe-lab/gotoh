#TODO:  Scoring algorithm to improve frame_align?

require 'cfe_gotoh/cfe_gotoh'


if(RUBY_PLATFORM =~ /(win|w)32$/)
  if(RUBY_VERSION =~ /^2/)
    require_relative 'alignment.windows.r22.so'
  else
    class String
      def ord() #adding an ord method
        return self.unpack('c')[0]
      end
    end
    require 'ckwlib/alignment.windows.so'
  end
elsif(RUBY_PLATFORM =~ /x86_64-linux/) #Ok, its probably not a mac, soo.....
  if(RUBY_VERSION =~ /^2/)
    require_relative 'alignment.linux64.r2.so'
  else
    class String
      def ord() #adding an ord method
        return self.unpack('c')[0]
      end
    end
    require 'ckwlib/alignment.linux64.so'
  end
elsif(RUBY_PLATFORM =~ /i686-darwin10/)
  require 'ckwlib/alignment.macosx.so'
else
  if(RUBY_VERSION =~ /^2/)
    require_relative 'alignment.linux32.r2.so'
  else
    class String
      def ord() #adding an ord method
        return self.unpack('c')[0]
      end
    end
    require 'ckwlib/alignment.linux32.so'    
  end
end

#init $nucMat
$nucMat = Array.new(127) {Array.new(127) {-1.0} }
['A','T','G','C','R','Y','K','M','B','D','H','V','S','W','N'].each do |nuc|
  $nucMat[nuc.ord()][nuc.ord()] = 1.0
  $nucMat[nuc.ord()]['X'.ord()]=$nucMat['X'.ord()][nuc.ord()]=-6.0 if(nuc !='N')
end
#bi-mixtures
$nucMat['A'.ord()]['R'.ord()]=$nucMat['R'.ord()]['A'.ord()]=1.0
$nucMat['G'.ord()]['R'.ord()]=$nucMat['R'.ord()]['G'.ord()]=1.0
$nucMat['C'.ord()]['Y'.ord()]=$nucMat['Y'.ord()]['C'.ord()]=1.0
$nucMat['T'.ord()]['Y'.ord()]=$nucMat['Y'.ord()]['T'.ord()]=1.0
$nucMat['G'.ord()]['K'.ord()]=$nucMat['K'.ord()]['G'.ord()]=1.0
$nucMat['T'.ord()]['K'.ord()]=$nucMat['K'.ord()]['T'.ord()]=1.0
$nucMat['C'.ord()]['M'.ord()]=$nucMat['M'.ord()]['C'.ord()]=1.0
$nucMat['A'.ord()]['M'.ord()]=$nucMat['M'.ord()]['A'.ord()]=1.0
$nucMat['C'.ord()]['S'.ord()]=$nucMat['S'.ord()]['C'.ord()]=1.0
$nucMat['G'.ord()]['S'.ord()]=$nucMat['S'.ord()]['G'.ord()]=1.0
$nucMat['T'.ord()]['W'.ord()]=$nucMat['W'.ord()]['T'.ord()]=1.0
$nucMat['A'.ord()]['W'.ord()]=$nucMat['W'.ord()]['A'.ord()]=1.0
#tri-mixtures
$nucMat['C'.ord()]['B'.ord()]=$nucMat['B'.ord()]['C'.ord()]=1.0
$nucMat['G'.ord()]['B'.ord()]=$nucMat['B'.ord()]['G'.ord()]=1.0
$nucMat['T'.ord()]['B'.ord()]=$nucMat['B'.ord()]['T'.ord()]=1.0
$nucMat['A'.ord()]['D'.ord()]=$nucMat['D'.ord()]['A'.ord()]=1.0
$nucMat['G'.ord()]['D'.ord()]=$nucMat['D'.ord()]['G'.ord()]=1.0
$nucMat['T'.ord()]['D'.ord()]=$nucMat['D'.ord()]['T'.ord()]=1.0
$nucMat['A'.ord()]['H'.ord()]=$nucMat['H'.ord()]['A'.ord()]=1.0
$nucMat['C'.ord()]['H'.ord()]=$nucMat['H'.ord()]['C'.ord()]=1.0
$nucMat['T'.ord()]['H'.ord()]=$nucMat['H'.ord()]['T'.ord()]=1.0
$nucMat['A'.ord()]['V'.ord()]=$nucMat['V'.ord()]['A'.ord()]=1.0
$nucMat['C'.ord()]['V'.ord()]=$nucMat['V'.ord()]['C'.ord()]=1.0
$nucMat['G'.ord()]['V'.ord()]=$nucMat['V'.ord()]['G'.ord()]=1.0
#other
$nucMat['$'.ord()]['$'.ord()]=50.0
$nucMat['T'.ord()]['U'.ord()] = $nucMat['U'.ord()]['T'.ord()] = 1.0
$nucMat['N'.ord()]['N'.ord()] = 0.0
$nucMat['X'.ord()]['-'.ord()]=$nucMat['X'.ord()]['-'.ord()]=3.0
['A','T','G','C'].each do |ch|
  $nucMat[ch.ord()]['*'.ord()]=$nucMat['*'.ord()][ch.ord()]=1.0
  $nucMat[ch.ord()]['&'.ord()]=$nucMat['&'.ord()][ch.ord()]=0.7
  $nucMat[ch.ord()]['$'.ord()]=$nucMat['$'.ord()][ch.ord()]=0.0
  $nucMat[ch.ord()]['.'.ord()]=$nucMat['.'.ord()][ch.ord()]=-20.0
  $nucMat[ch.ord()]['N'.ord()]=$nucMat['N'.ord()][ch.ord()]=-3.0
end



def score_alignment(seqa, seqb)
  sc = 0.0
  0.upto(seqa.size() - 1) do |i|
    sc += $nucMat[seqa[i,1].upcase().ord()][seqb[i,1].upcase().ord()]
  end
  return sc
end

def make_gap_list(seq)
  list = []
  cur_ins = nil
  prev_i = nil
  0.upto(seq.size() - 1) do |i|
    if(seq[i,1] == '-')
      if(prev_i and i == prev_i + 1)
        cur_ins << i
        prev_i = i
      else
        list << cur_ins if(cur_ins != nil and cur_ins != [])
        cur_ins = [i]
        prev_i = i
      end
    end
  end
  list << cur_ins if(cur_ins != nil and cur_ins != [])
  return list
end

#common_insert_locations is based on amino acid locations starting at base 0.
#Assumes standard in the first base.
#Prealign lets you run a lot of the corrections and qc on a already aligned sequence.  
def frame_align(seqa, seqb, gap_init=3, gap_penalty=1, common_insert_locations=[], trim=false, raise_errors=false, prealigned=false)
  elem = nil
  if(prealigned)
    elem = [seqa, seqb]
  else
    elem = align_it(seqa, seqb, gap_init, gap_penalty)
  end
  puts "Wierd sizes Z" if(elem[0].size() != elem[1].size())
  
  #Do trimming?
  if(trim and elem[0] =~ /^(-+)[^-]/)
    elem[0][0,$1.size()] = ''
    elem[1][0,$1.size()] = ''
  end
  if(trim and elem[0] =~ /[^-](-+)$/)
    elem[1][(elem[0].size() - $1.size()), $1.size()] = ''
    elem[0][(elem[0].size() - $1.size()), $1.size()] = ''
  end
  
  #Start 
  if(trim and elem[1][0,1] == '-')
    #get rid of edges that are the wrong size
    elem[1] =~ /^(-+)[^-]/
    #Make sure its a multiple of three
    dashes = $1
    if(dashes == nil)
      #pass
    elsif((dashes.size() % 3) == 1) 
      elem[1][dashes.size(),1] = '-'
      elem[1][dashes.size() + 1,1] = '-'
    elsif((dashes.size() % 3) == 2)
      elem[1][dashes.size(),1] = '-'
    end  
  end
  
  #end
  if(trim and elem[1][-1,1] == '-')
    #get rid of edges that are the wrong size
    elem[1] =~ /[^-](-+)$/
    #Make sure its a multiple of three
    dashes = $1
    if(dashes == nil)
      #pass
    elsif((dashes.size() % 3) == 1)
      elem[1][(elem[1].size() - dashes.size()) - 1] = '-'
      elem[1][(elem[1].size() - dashes.size()) - 2] = '-'
    elsif((dashes.size() % 3) == 2)
      elem[1][(elem[1].size() - dashes.size()) - 1] = '-'
    end  
  end
  
  
  #try to merge deletions and insertions if things aren't looking well.
  #added 16-Nov-2018, helps fix poor insertions near the start.
  if(elem[0].size() % 3 != 0 or elem[1].size() % 3 != 0)
    dex = 0 #Don't start at 0, that way lies madness...  Or does it???
    while(dex = elem[0].index(/-/, dex))
      #Now, find an ajacent dash in elem[1] to cancel
      if((dex - 1 >= 0) and elem[1][dex - 1] == '-')
        elem[0][dex] = ''
        elem[1][dex - 1] = ''
        dex = 1
      elsif(elem[1][dex + 1] == '-')
        elem[1][dex + 1] = ''
        elem[0][dex] = ''
        dex = 1
      elsif((dex - 2 >= 0) and elem[1][dex - 2] == '-')
        elem[1][dex - 2] = ''
        elem[0][dex] = ''
        dex = 1
      elsif(elem[1][dex + 2] == '-')
        elem[1][dex + 2] = ''
        elem[0][dex] = ''
        dex = 1
      end
      
      #check to see if we fixed everything
      if(!(elem[0].size() % 3 != 0 or elem[1].size() % 3 != 0))
        break
      end
    
      dex += 1
    end
    
  end

  
  #I wonder if these should throw exceptions by default, but have an option to ignore.
  if(elem[0].gsub(/[^-]/,'').size() % 3 != 0 and raise_errors)
    #puts "Can not frame align"
    raise "Can not frame align, #{elem[0].gsub(/[^-]/,'').size()} inserted bases not divisible by 3"
    #return elem
  end
  if(elem[1].gsub(/[^-]/,'').size() % 3 != 0 and raise_errors)
    #puts "Can not frame align"
    #puts elem[0]
    #puts elem[1]
    
    raise "Can not frame align, #{elem[1].gsub(/[^-]/,'').size()} deleted bases not divisible by 3"
    #return elem
  end
  
  #Build the insert/delete lists.
  insert_list = make_gap_list(elem[0])
  delete_list = make_gap_list(elem[1])
  #Now we have a list that looks like [[3,4,5], [9], [11,12]]
  

  if(insert_list.size() > 0)#Inserts first
    new_ins_list = []
    
    #First step is clustering insertions. (v2:  16-Nov-2018)
    insert_list.each_with_index do |ins, i|
      next if(ins.size() == 0) #we already ate this one.
      if(ins.size() % 3 == 0) #this insertion is fine!
        new_ins_list << ins
        next 
      end
      
      #Can I merge with the next insert?
      if(insert_list[i + 1] and (ins + insert_list[i + 1]).size() % 3 == 0 and
        (insert_list[i + 1].first - ins.last) < 9)  
        
        ins2 = insert_list[i + 1]
        if(ins2.size() > ins.size())
          new_ins_list <<  ((ins2.first - ins.size()) .. ins2.first - 1).to_a() + ins2
        else
          new_ins_list <<  ins + ((ins.last + 1) .. (ins.last + ins2.size())).to_a()
        end
        insert_list[i + 1] = []
      #maybe merge with the next two inserts?
      elsif(insert_list[i + 1] and insert_list[i + 2] and
        (ins + insert_list[i + 1] + insert_list[i + 2]).size() % 3 == 0 and
        (insert_list[i + 2].first - ins.last) < 12) 
        
        ins2 = insert_list[i + 1]
        ins3 = insert_list[i + 2]
        if(true) #Lets just assume that if you need to combine 3 inserts, the middle one is where it goes.
          new_ins_list << ((ins2.first - ins.size()) .. ins2.first - 1).to_a() + ins2 + ((ins2.last + 1) .. (ins2.last + ins3.size())).to_a()
        end
        
        insert_list[i + 1] = []
        insert_list[i + 2] = []
      else #No merge, life sucks and then you die.
        raise "Can not frame align insert" if(raise_errors)
      end
    end
    
=begin
    #First step is clustering insertions. (v1:  old version)
    insert_list.each_with_index do |ins, i|
      next_ins = insert_list[i + 1]
      
      if(ins.size() % 3 != 0) #Wrong size!
        #Look for next insertions that would make it the right size and are not far apart
        if(!outta_frame and next_ins and (((next_ins.size() + ins.size()) % 3) == 0) and next_ins[0] - ins[-1] < 8 ) #within 8 bases
          if(next_ins.size() > ins.size()) #scoring would be good here
            ins.each do |a|
              next_ins.insert(0, next_ins[0] - 1) #Insert at start
            end
          else
            next_ins.each do |a|
              ins.insert(-1, ins[0] + 1) #insert at end
            end
            insert_list[i + 1] = ins #Pushing this problem ahead for convinence.
          end
          #don't insert into new list, as we've pushed it to the next element.
        else
          #I wonder if we should try to merge multiple more than two inserts?  Eh, nah.
          #return elem
          outta_frame = true
          raise "Can not frame align insert" if(raise_errors)
        end
      else
        new_ins_list << ins #this insert is okay
      end
      
    end
=end
    
    #puts insert_list.inspect
    #puts new_ins_list.inspect
    
    #second step should be to frame align inserts (prioritizing to common_points)
    offset = 0 #offset created by previous insertions. (IMPORTANT)
    new_ins_list.each do |ins|
      #see if its close to a common_insert(within 3 amino acids?)
      min_common = common_insert_locations.min(){|a,b| ((a) * 3 - (ins[0] - offset)).abs() <=> ((b) * 3 - ins[0]).abs()}
      if(min_common != nil and ((min_common ) * 3 - (ins[0] - offset)).abs() <= 9)
        #Cool, align to this common insert
        new_ins = []
        0.upto(ins.size() - 1) do |i|
          new_ins << ((min_common) * 3) + i + offset
        end
        ins.replace(new_ins)
      end
      
      #B frame align 
      #scoring would be good here
      if(ins[0] % 3 == 1) #set back one base
        new_ins = []
        ins.each do |i|
          new_ins << i - 1
        end
        ins.replace(new_ins)
      elsif(ins[0] % 3 == 2) #Set forward one base.
        new_ins = []
        ins.each do |i|
          new_ins << i + 1
        end
        ins.replace(new_ins)
      end
      
      offset += ins.size()
    end

    #make the actual modifications
    #begin
    #orige = "" + elem[0]
    elem[0] = elem[0].gsub('-','')
    new_ins_list.each do |ins|
      ins.each do |i|
        if(i > elem[0].size())
          elem[0].insert(-1, '-')
        else
          elem[0].insert(i, '-')
        end
      end
    end
    #rescue
    #  puts "OOH---------------------------------------------"
    #  puts orige.inspect
    #  puts elem.inspect
    #  puts new_ins_list.inspect()
    #  raise $!
    #end

  end
  
  
  #Deletion---------------------------------------------------------------------------------------------------
  outta_frame = false
  if(delete_list.size() > 0)#Deletions second
    new_del_list = []
    next_del = nil
    
    #First step is clustering deletions. (v2:  19-Nov-2018)
    delete_list.each_with_index do |del, i|
      next if(del.size() == 0) #we already ate this one.
      if(del.size() % 3 == 0) #this insertion is fine!
        new_del_list << del
        next 
      end
      
      #Can I merge with the next delete?
      if(delete_list[i + 1] and (del + delete_list[i + 1]).size() % 3 == 0 and
        (delete_list[i + 1].first - del.last) < 9)  
        
        del2 = delete_list[i + 1]
        if(del2.size() > del.size())
          new_del_list <<  ((del2.first - del.size()) .. del2.first - 1).to_a() + del2
        else
          new_del_list <<  del + ((del.last + 1) .. (del.last + del2.size())).to_a()
        end
        delete_list[i + 1] = []
      #maybe merge with the next two deletes?
      elsif(delete_list[i + 1] and delete_list[i + 2] and
        (del + delete_list[i + 1] + delete_list[i + 2]).size() % 3 == 0 and
        (delete_list[i + 2].first - del.last) < 12)  #slightly higher range, since we've already got a higher threshold of confidence here
        
        del2 = delete_list[i + 1]
        del3 = delete_list[i + 2]
        if(true) #Lets just assume that if you need to combine 3 deletes, the middle one is where it goes.
          new_del_list << ((del2.first - del.size()) .. del2.first - 1).to_a() + del2 + ((del2.last + 1) .. (del2.last + del3.size())).to_a()
        end
        
        delete_list[i + 1] = []
        delete_list[i + 2] = []
      else #No merge, life sucks and then you die.
        new_del_list << delete_list[i]
      end
    end
    
=begin
    #First step is clustering deletions. (v1:  old version)
    delete_list.each_with_index do |del, i|
      next_del = delete_list[i + 1]
      
      if(del.size() % 3 != 0 and !outta_frame) #Wrong size!
        #Look for next deletions that would make it the right size and are not far apart
        if(next_del and (((next_del.size() + del.size()) % 3) == 0) and next_del[0] - del[-1] < 6 ) #within 6 bases
          if(next_del.size() > del.size()) #scoring would be good here
            del.each do |a|
              next_del.insert(0, next_del[0] - 1) #delete at start
            end
          else
            next_del.each do |a|
              del.insert(-1, del[0] + 1) #delete at end
            end
            delete_list[i + 1] = del #Pushing this problem ahead for convinence.
          end
          #don't delete into new list, as we've pushed it to the next element.
        else
          #I wonder if we should try to merge multiple more than two deletes?  Eh, nah.
          #
          #return elem
          new_del_list << del 
          outta_frame = true
          #raise "Can not frame align delete"
        end
      else
        new_del_list << del #this delete is okay
      end
    end
=end

    
    #second step should be to frame align deletes 
    offset = 0 #offset created by previous deletions. (IMPORTANT)
    new_del_list.each do |del|
      next if(del.size() % 3 != 0)
      #frame align 
      #scoring would be good here
      if(del[0] % 3 == 1) #set back one base
        new_del = []
        del.each do |i|
          new_del << i - 1
        end
        del.replace(new_del)
      elsif(del[0] % 3 == 2) #Set forward one base.
        new_del = []
        del.each do |i|
          new_del << i + 1
        end
        del.replace(new_del)
      end
      
      offset += del.size()
    end
    
    #make the actual modifications
    elem[1] = elem[1].gsub('-','')
    new_del_list.each do |del|
      del.each do |i|
        if(i > elem[1].size() )
          elem[1].insert(elem[1].size(), '-')
        else
          elem[1].insert(i, '-')
        end
      end
    end
  end

  return elem
end

#Returns a [seq_sans_inserts, [list of inserts]]
def remove_inserts(elem)
  seq = '' + elem[1]
  inserts = []
  
  insert_list = []
  0.upto(elem[0].size() - 1) do |i|
    insert_list << i if(elem[0][i,1] == '-')
  end
  
  big_insert_list = []
  if(elem[0].include?('-'))#Inserts first
    #First step should be to cluster inserts
    cur_ins = nil
    prev_i = nil
    insert_list.each do |i|
      if(prev_i and i == prev_i + 1)
        cur_ins << i
        prev_i = i
      else
        big_insert_list << cur_ins if(cur_ins != nil and cur_ins != [])
        cur_ins = [i]
        prev_i = i
      end
    end
    big_insert_list << cur_ins if(cur_ins != nil and cur_ins != [])
  end
  
  offset = 0
  big_insert_list.each do |ins|
    ins_seq = ''
    ins.each do |i|
      ins_seq += elem[1][i,1]
    end
    inserts << [((ins[0] - offset) / 3), ins_seq]
    offset += ins.size()
    ins.each do |i|
      seq[i,1] = '.'
    end
  end

  return [seq.gsub('.',''), inserts]
end

=begin
seq = 'CCTCAAATCACTCTTTGGCAACGACCCTTAGTCACAGTAAGAATAGGGGGACAGCTAATAGAAGCCCTATTAGACACAGGAGCAGATGATACAGTATTAGAAGAAAAAATAGATTTACCAGGAAAATGGARACCAAAAATGATAGGGGGAATTGGAGGTTTTATTAAAGTAAGGCAATATGATCAGATACTTATGGAAATATGTGARAAGAAGGCCATAGGTACAGTATTAGTAGGACCTACMCCTGTCAACATAATTGGRCGRAATATGTTGACTCAGATTGGTTGTACTTTAAATTTTCCAATTAGTCCTATTGARACTGTGCCAGTAAAATTAAAGCCAGGGATGGATGGCCCAAAAGTTAARCAATGGCCATTGACAGAAGAAAAAATAAAAGCATTAACAGAAATTTGTGCAGAAATGGAAAARGAAGGAAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTGTTTGCCATAAAGAAAAARGATAGTAMTAAATGGAGAAAATTAGTAGATTTCAGAGAACTCAATAAGAGAACTCAAGACTTCTGGGAGGTCCAATTAGGAATTCCTCATCCCGCGGGATTAAAAAAGAAAARATCAGTAACAGTACTAGATRTAGGGGATGCATATTTTTCAGTTCCCTTAGACAAAGAYTTTAGAAAGTATACTGCATTCACTATACCTAGTGTAAATAATGAAACACCAGGRATTAGATATCAGTACAATGTRCTKCCACAGGGATGGAAAGGATCACCAGCAATATTTCARGCAAGCATGACAAAAATCTTAGAGCCCTTTAGAACAAAAAATCCAGAGGTGGTGATCTACCAGTATATGGATGATTTATATGTAGGATCTGACTTAGAGATAGGGCAACATAGAGCAAAAATAGAGGARTTAAGAGAACATCTAYTGARATGGGGATTTACCACACCAGACAAAAAACATCAGAAAGAACCTCCATTTCTTTGGATGGGATATGAACTTCATCCTGACAAATGGACAGTCCAGCCTATARTRCTGCCARACAAAGRMRRCTGGACTGTCAATGATATACAGAAATTAGTAGGAAAACTAAATTGGGCCAGTCAAATTTATGCAGGAATTAAAGTAAAGCAACTGTGTAAACTCCTCAGGGGAGCCAAAGCATTAACAGAYATAGTAACAYTAACTGAGGAAGCAGAATTAGAATTGGCAGAGAACAGGGAAATTCTAAAAGAACCTGTACATGGGGTATAYTATGAYCCAGYAAAAGACTTAATAGCAGAAATACAGAAACAAGGGCAAGACCAATGGACATATCAAATATATCAAGARCCATTTAAAAATCTAAARACAGGAAAATATGCAAARAGGAGATCTGCCCACACRAATGATGTAAAACAATTAACAGAGGTAGTGCAAAAAGTGTCTACAGAARGCATAGTAATATGGGGGAAGAYCCCTAAATTTAAGCTGCCCATACAAAAAGAAACATGGGAGGCA'
std = 'CCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTTCCCATTAGCCCTATTGAGACTGTACCAGTAAAATTAAAGCCAGGAATGGATGGCCCAAAAGTTAAACAATGGCCATTGACAGAAGAAAAAATAAAAGCATTAGTAGAAATTTGTACAGAGATGGAAAAGGAAGGGAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAAGAAAAAAGACAGTACTAAATGGAGAAAATTAGTAGATTTCAGAGAACTTAATAAGAGAACTCAAGACTTCTGGGAAGTTCAATTAGGAATACCACATCCCGCAGGGTTAAAAAAGAAAAAATCAGTAACAGTACTGGATGTGGGTGATGCATATTTTTCAGTTCCCTTAGATGAAGACTTCAGGAAGTATACTGCATTTACCATACCTAGTATAAACAATGAGACACCAGGGATTAGATATCAGTACAATGTGCTTCCACAGGGATGGAAAGGATCACCAGCAATATTCCAAAGTAGCATGACAAAAATCTTAGAGCCTTTTAGAAAACAAAATCCAGACATAGTTATCTATCAATACATGGATGATTTGTATGTAGGATCTGACTTAGAAATAGGGCAGCATAGAACAAAAATAGAGGAGCTGAGACAACATCTGTTGAGGTGGGGACTTACCACACCAGACAAAAAACATCAGAAAGAACCTCCATTCCTTTGGATGGGTTATGAACTCCATCCTGATAAATGGACAGTACAGCCTATAGTGCTGCCAGAAAAAGACAGCTGGACTGTCAATGACATACAGAAGTTAGTGGGGAAATTGAATTGGGCAAGTCAGATTTACCCAGGGATTAAAGTAAGGCAATTATGTAAACTCCTTAGAGGAACCAAAGCACTAACAGAAGTAATACCACTAACAGAAGAAGCAGAGCTAGAACTGGCAGAAAACAGAGAGATTCTAAAAGAACCAGTACATGGAGTGTATTATGACCCATCAAAAGACTTAATAGCAGAAATACAGAAGCAGGGGCAAGGCCAATGGACATATCAAATTTATCAAGAGCCATTTAAAAATCTGAAAACAGGAAAATATGCAAGAATGAGGGGTGCCCACACTAATGATGTAAAACAATTAACAGAGGCAGTGCAAAAAATAACCACAGAAAGCATAGTAATATGGGGAAAGACTCCTAAATTTAAACTGCCCATACAAAAGGAAACATGGGAAACATGGTGGACAGAGTATTGGCAAGCCACCTGGATTCCTGAGTGGGAGTTTGTTAATACCCCTCCCTTAGTGAAATTATGGTACCAGTTAGAGAAAGAACCCATAGTAGGAGCAGAAACCTTC'

elem = frame_align(std, seq, 9, 2, [35,168]) #Seems to be off by one.
puts remove_inserts(elem).inspect()
=end
=begin
#----testing-----
require 'ckwlib/io'

std = "TGTACAAGACCCAACAACAATACAAGAAAAAGTATACATATAGGACCAGGGAGAGCATTTTATGCAACAGGAGAAATAATAGGAGATATAAGACAAGCACATTGT"

fasta = Io::read_fasta('motivate_preds_full_NucSeq.fas')

File.open(

fasta.each do |fas|
  elem = frame_align(std, fas[1].strip(), 3, 1, [10])
  puts fas[0]
  puts elem[0]
  puts elem[1]
  puts remove_inserts(elem).inspect()
  puts score_alignment(elem[0], elem[1])
end
=end