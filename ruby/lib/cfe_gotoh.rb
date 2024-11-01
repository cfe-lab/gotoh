#TODO:  Scoring algorithm to improve frame_align?

require 'cfe_gotoh/cfe_gotoh'


module CfeGotoh
  class Error < RuntimeError
  end

  class GapMergeError < Error
  end

  sub_matrix = Array.new(127) {Array.new(127) {-1.0} }
  ['A','T','G','C','R','Y','K','M','B','D','H','V','S','W','N'].each do |nuc|
    sub_matrix[nuc.ord()][nuc.ord()] = 1.0
    sub_matrix[nuc.ord()]['X'.ord()]=sub_matrix['X'.ord()][nuc.ord()]=-6.0 if(nuc !='N')
  end
  #bi-mixtures
  sub_matrix['A'.ord()]['R'.ord()]=sub_matrix['R'.ord()]['A'.ord()]=1.0
  sub_matrix['G'.ord()]['R'.ord()]=sub_matrix['R'.ord()]['G'.ord()]=1.0
  sub_matrix['C'.ord()]['Y'.ord()]=sub_matrix['Y'.ord()]['C'.ord()]=1.0
  sub_matrix['T'.ord()]['Y'.ord()]=sub_matrix['Y'.ord()]['T'.ord()]=1.0
  sub_matrix['G'.ord()]['K'.ord()]=sub_matrix['K'.ord()]['G'.ord()]=1.0
  sub_matrix['T'.ord()]['K'.ord()]=sub_matrix['K'.ord()]['T'.ord()]=1.0
  sub_matrix['C'.ord()]['M'.ord()]=sub_matrix['M'.ord()]['C'.ord()]=1.0
  sub_matrix['A'.ord()]['M'.ord()]=sub_matrix['M'.ord()]['A'.ord()]=1.0
  sub_matrix['C'.ord()]['S'.ord()]=sub_matrix['S'.ord()]['C'.ord()]=1.0
  sub_matrix['G'.ord()]['S'.ord()]=sub_matrix['S'.ord()]['G'.ord()]=1.0
  sub_matrix['T'.ord()]['W'.ord()]=sub_matrix['W'.ord()]['T'.ord()]=1.0
  sub_matrix['A'.ord()]['W'.ord()]=sub_matrix['W'.ord()]['A'.ord()]=1.0
  #tri-mixtures
  sub_matrix['C'.ord()]['B'.ord()]=sub_matrix['B'.ord()]['C'.ord()]=1.0
  sub_matrix['G'.ord()]['B'.ord()]=sub_matrix['B'.ord()]['G'.ord()]=1.0
  sub_matrix['T'.ord()]['B'.ord()]=sub_matrix['B'.ord()]['T'.ord()]=1.0
  sub_matrix['A'.ord()]['D'.ord()]=sub_matrix['D'.ord()]['A'.ord()]=1.0
  sub_matrix['G'.ord()]['D'.ord()]=sub_matrix['D'.ord()]['G'.ord()]=1.0
  sub_matrix['T'.ord()]['D'.ord()]=sub_matrix['D'.ord()]['T'.ord()]=1.0
  sub_matrix['A'.ord()]['H'.ord()]=sub_matrix['H'.ord()]['A'.ord()]=1.0
  sub_matrix['C'.ord()]['H'.ord()]=sub_matrix['H'.ord()]['C'.ord()]=1.0
  sub_matrix['T'.ord()]['H'.ord()]=sub_matrix['H'.ord()]['T'.ord()]=1.0
  sub_matrix['A'.ord()]['V'.ord()]=sub_matrix['V'.ord()]['A'.ord()]=1.0
  sub_matrix['C'.ord()]['V'.ord()]=sub_matrix['V'.ord()]['C'.ord()]=1.0
  sub_matrix['G'.ord()]['V'.ord()]=sub_matrix['V'.ord()]['G'.ord()]=1.0
  #other
  sub_matrix['$'.ord()]['$'.ord()]=50.0
  sub_matrix['T'.ord()]['U'.ord()] = sub_matrix['U'.ord()]['T'.ord()] = 1.0
  sub_matrix['N'.ord()]['N'.ord()] = 0.0
  sub_matrix['X'.ord()]['-'.ord()]=sub_matrix['X'.ord()]['-'.ord()]=3.0
  ['A','T','G','C'].each do |ch|
    sub_matrix[ch.ord()]['*'.ord()]=sub_matrix['*'.ord()][ch.ord()]=1.0
    sub_matrix[ch.ord()]['&'.ord()]=sub_matrix['&'.ord()][ch.ord()]=0.7
    sub_matrix[ch.ord()]['$'.ord()]=sub_matrix['$'.ord()][ch.ord()]=0.0
    sub_matrix[ch.ord()]['.'.ord()]=sub_matrix['.'.ord()][ch.ord()]=-20.0
    sub_matrix[ch.ord()]['N'.ord()]=sub_matrix['N'.ord()][ch.ord()]=-3.0
  end
  sub_matrix.each {|column| column.freeze}
  sub_matrix.freeze

  NUCLEOTIDE_MATRIX = sub_matrix

  def self.score_alignment(standard, query)
    sc = 0.0
    0.upto(standard.size() - 1) do |i|
      sc += NUCLEOTIDE_MATRIX[standard[i,1].upcase().ord()][query[i,1].upcase().ord()]
    end
    return sc
  end

  def self.make_gap_list(seq)
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

  def self.trim_leading_dashes(standard, query)
    leading_dashes_match = /^(-+)[^-]/.match(standard)
    if (leading_dashes_match.nil?)
      return
    end
    leading_dashes = leading_dashes_match[1]
    standard[0, leading_dashes.size()] = ''
    query[0, leading_dashes.size()] = ''
  end

  def self.trim_trailing_dashes(standard, query)
    trailing_dashes_match = /[^-](-+)$/.match(standard)
    if (trailing_dashes_match.nil?)
      return
    end
    trailing_dashes = trailing_dashes_match[1]
    end_of_standard = standard.size() - trailing_dashes.size()
    standard[end_of_standard, trailing_dashes.size()] = ''
    query[end_of_standard, trailing_dashes.size()] = ''
  end

  def self.fix_incomplete_edge_codon(query, side=:leading)
    edge_idx = 0
    dash_regex = /^(-+)[^-]/
    incr = 1
    if (side != :leading)  # fix the trailing edge
      edge_idx = -1
      dash_regex = /[^-](-+)$/
      incr = -1
    end

    if (query[edge_idx] == '-')
      dashes = dash_regex.match(query)[0]  # we know there will be a match

      # If the length of the dashes aren't a multiple of 3, turn some
      # of the query characters into dashes to force it to be a full
      # codon of dashes.
      if (dashes.size() % 3 >= 1)
        first_non_dash_idx = 0
        if (side != :leading)
          first_non_dash_idx = query.size() - dashes.size() - 1
        end
        query[first_non_dash_idx] = '-'
        if (dashes.size() % 3 == 2)
          query[first_non_dash_idx + incr] = '-'
        end
      end
    end
  end

  def self.merge_insertions_and_deletions_to_fix_oof_sequences(
    standard,
    query
  )
    # Merge deletions and insertions until the sequences have a cogent length
    # (i.e. have length divisible by 3).  This helps fix poor insertions near 
    # the start of the sequence.
    raise 'Standard and query should be the same length' if standard.size() != query.size()
    if(standard.size() % 3 != 0)
      dex = 0
      while(dex = standard.index(/-/, dex))
        [-1, 1, -2, 2].each do |offset|  # look one base away, then two bases away
          if ((dex + offset >= 0) and query[dex + offset] == '-')
            standard[dex] = ''
            query[dex + offset] = ''
            dex = 0
            break
          end
        end
        
        # Stop if the sequences are now a cogent length.
        if(standard.size() % 3 == 0)
          break
        end
        dex += 1
      end
    end
  end

  def self.cluster_gaps(gaps, raise_errors=false)
    # Merge adjacent gaps if they are not a codon-sized gap.
    new_gap_list = []
    gaps.each_with_index do |gap, i|
      next if(gap.size() == 0)  # we already ate this one
      if(gap.size() % 3 == 0)  # this gap is fine!
        new_gap_list << gap
        next
      end
      
      gap2 = gaps[i + 1]  # note: these could be nil, which we test for below
      gap3 = gaps[i + 2]
      # Can I merge with the next gap?
      if (gap2 and (gap + gap2).size() % 3 == 0 and (gap2.first - gap.last) < 9)
        if(gap2.size() > gap.size())
          new_gap_list <<  ((gap2.first - gap.size()) .. gap2.first - 1).to_a() + gap2
        else
          new_gap_list <<  gap + ((gap.last + 1) .. (gap.last + gap2.size())).to_a()
        end
        gaps[i + 1] = []
      # Can I merge with the next two gaps?
      elsif(
        gap2 and gap3 and
        (gap + gap2 + gap3).size() % 3 == 0 and
        (gap3.first - gap.last) < 12
      )
        # Place the gap around the middle of the three merging gaps.
        new_gap = (
          (gap2.first - gap.size()) .. gap2.first - 1.to_a()
          + gap2
          + ((gap2.last + 1) .. (gap2.last + gap3.size())).to_a()
        )
        new_gap_list << new_gap
        
        gaps[i + 1] = []
        gaps[i + 2] = []
      else
        # We can't merge the gaps; either raise an error or meekly proceed.
        if (raise_errors)
          raise GapMergeError
        else
          new_gap_list << gap  # FIXME this behaviour differs between insertions and deletions
        end
      end
    end
    return new_gap_list
  end

  def self.align_gaps_to_frame(gaps, common_gap_locations=nil)
    # Align gaps to codon boundaries, giving preference to common 
    # gap locations if specified.
    
    offset = 0  # offset created by previous gaps.
    gaps.each do |gap|
      # See if this gap is close to a common gap location (within 3 amino acids).
      if (!common_gap_locations.nil?)
        closest_common = common_gap_locations.min() do |a,b|
          (3 * a - (gap[0] - offset)).abs() <=> (3 * b - (gap[0] - offset)).abs()
        end
        if(closest_common != nil and (3 * closest_common - (gap[0] - offset)).abs() <= 9)
          # Align the gap to this position.
          new_gap = []
          0.upto(gap.size() - 1) do |i|
            new_gap << 3 * closest_common + i + offset
          end
          gap.replace(new_gap)
        end
      end
      
      # Align the gap to the nearest appropriate frame.
      # Original comment from Conan: scoring would be good here
      if(gap[0] % 3 == 1)  # set back one base
        new_gap = []
        gap.each do |i|
          new_gap << i - 1
        end
        gap.replace(new_gap)
      elsif(gap[0] % 3 == 2)  # set forward one base
        new_gap = []
        gap.each do |i|
          new_gap << i + 1
        end
        gap.replace(new_gap)
      end
      
      offset += gap.size()
    end
    return gaps
  end

  def self.splice_gaps_into_sequence(seq, gaps)
    # Place the specified gaps into the sequence.  Note that the
    # gaps are specified by their positions in an *aligned* sequence,
    # and as such include "offsets" introduced by gaps placed earlier
    # in the sequence.  The gaps must be in left-to-right order.
    seq = seq.gsub('-','')
    gaps.each do |gap|
      gap.each do |i|
        if(i > seq.size())
          seq.insert(-1, '-')
        else
          seq.insert(i, '-')
        end
      end
    end
    return seq
  end

  #common_insert_locations is based on amino acid locations starting at base 0.
  #Assumes standard in the first base.
  #Prealign lets you run a lot of the corrections and qc on a already aligned sequence.  
  def self.frame_align(
    standard,
    query,
    gap_init=3,
    gap_penalty=1,
    common_insert_locations=nil,
    trim=false,
    raise_errors=false,
    prealigned=false
  )
    if (common_insert_locations.nil?)
      common_insert_locations = []
    end
    if(!prealigned)
      elem = align_it(standard, query, gap_init, gap_penalty)
      standard = elem[0]
      query = elem[1]
    end
    raise "Standard and query should be the same length" if standard.size() != query.size()
    
    # Trim leading and trailing dashes if desired.
    if (trim)
      trim_leading_dashes(standard, query)
      trim_trailing_dashes(standard, query)
      fix_incomplete_edge_codon(standard, query, :leading)
      fix_incomplete_edge_codon(standard, query, :trailing)
    end
    
    merge_insertions_and_deletions_to_fix_oof_sequences(standard, query)

    if(standard.gsub(/[^-]/,'').size() % 3 != 0 and raise_errors)
      raise "Cannot frame align, #{standard.gsub(/[^-]/,'').size()} inserted bases not divisible by 3"
    end
    if(query.gsub(/[^-]/,'').size() % 3 != 0 and raise_errors)
      raise "Cannot frame align, #{query.gsub(/[^-]/,'').size()} deleted bases not divisible by 3"
    end
    
    # Build the insert/delete lists.  These lists look like 
    # [[3,4,5], [9], [11,12]]
    insert_list = make_gap_list(standard)
    delete_list = make_gap_list(query)
    
    # Process the insertions.
    if(insert_list.size() > 0)
      new_ins_list = []
      
      # Step 1: cluster the insertions.
      begin
        new_ins_list = cluster_gaps(insert_list, raise_errors=raise_errors)
      rescue GapMergeError
        raise "Cannot frame align insert" if raise_errors
      end
      
      # Step 2: frame-align the insertions, shifting things to common insertion
      # positions where appropriate.
      align_gaps_to_frame(new_ins_list, common_gap_locations=common_insert_locations)

      # Put the insertions back into the standard.
      standard = splice_gaps_into_sequence(standard, new_ins_list)
    end
    
    # Process the deletions.
    if(delete_list.size() > 0)
      new_del_list = []
      
      # As above, step 1 is to cluster the deletions.
      # FIXME note that the original code behaved differently between
      # insertions and deletions; confirm that this is the right
      # way forward.
      begin
        new_del_list = cluster_gaps(delete_list, raise_errors=raise_errors)
      rescue GapMergeError
        raise "Cannot frame align deletion" if raise_errors
      end
      
      # Again as above, frame-align the deletions; this time
      # we don't worry about any common deletion positions.
      align_gaps_to_frame(new_del_list)

      # Put the deletions back into the query.
      query = splice_gaps_into_sequence(query, new_del_list)
    end

    return [standard, query]
  end

  #Returns a [seq_sans_inserts, [list of inserts]]
  def self.remove_inserts(elem)
    return remove_insertions_from_query(elem[1])
  end

  def self.remove_insertions_from_query(query)
    seq = '' + query
    inserts = []
    
    insert_list = []
    0.upto(standard.size() - 1) do |i|
      insert_list << i if(standard[i,1] == '-')
    end
    
    big_insert_list = []
    if(standard.include?('-'))#Inserts first
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
        ins_seq += query[i,1]
      end
      inserts << [((ins[0] - offset) / 3), ins_seq]
      offset += ins.size()
      ins.each do |i|
        seq[i,1] = '.'
      end
    end

    return [seq.gsub('.',''), inserts]
  end
end
