#TODO:  Scoring algorithm to improve frame_align?

require_relative 'cfe_gotoh/cfe_gotoh'


module CfeGotoh
  class Error < RuntimeError
  end

  class GapMergeError < Error
  end

  sub_matrix = Array.new(127) {Array.new(127) {-1.0} }
  ['A','T','G','C','R','Y','K','M','B','D','H','V','S','W','N'].each do |nuc|
    sub_matrix[nuc.ord][nuc.ord] = 1.0
    sub_matrix[nuc.ord]['X'.ord]=sub_matrix['X'.ord][nuc.ord]=-6.0 if(nuc !='N')
  end
  #bi-mixtures
  sub_matrix['A'.ord]['R'.ord]=sub_matrix['R'.ord]['A'.ord]=1.0
  sub_matrix['G'.ord]['R'.ord]=sub_matrix['R'.ord]['G'.ord]=1.0
  sub_matrix['C'.ord]['Y'.ord]=sub_matrix['Y'.ord]['C'.ord]=1.0
  sub_matrix['T'.ord]['Y'.ord]=sub_matrix['Y'.ord]['T'.ord]=1.0
  sub_matrix['G'.ord]['K'.ord]=sub_matrix['K'.ord]['G'.ord]=1.0
  sub_matrix['T'.ord]['K'.ord]=sub_matrix['K'.ord]['T'.ord]=1.0
  sub_matrix['C'.ord]['M'.ord]=sub_matrix['M'.ord]['C'.ord]=1.0
  sub_matrix['A'.ord]['M'.ord]=sub_matrix['M'.ord]['A'.ord]=1.0
  sub_matrix['C'.ord]['S'.ord]=sub_matrix['S'.ord]['C'.ord]=1.0
  sub_matrix['G'.ord]['S'.ord]=sub_matrix['S'.ord]['G'.ord]=1.0
  sub_matrix['T'.ord]['W'.ord]=sub_matrix['W'.ord]['T'.ord]=1.0
  sub_matrix['A'.ord]['W'.ord]=sub_matrix['W'.ord]['A'.ord]=1.0
  #tri-mixtures
  sub_matrix['C'.ord]['B'.ord]=sub_matrix['B'.ord]['C'.ord]=1.0
  sub_matrix['G'.ord]['B'.ord]=sub_matrix['B'.ord]['G'.ord]=1.0
  sub_matrix['T'.ord]['B'.ord]=sub_matrix['B'.ord]['T'.ord]=1.0
  sub_matrix['A'.ord]['D'.ord]=sub_matrix['D'.ord]['A'.ord]=1.0
  sub_matrix['G'.ord]['D'.ord]=sub_matrix['D'.ord]['G'.ord]=1.0
  sub_matrix['T'.ord]['D'.ord]=sub_matrix['D'.ord]['T'.ord]=1.0
  sub_matrix['A'.ord]['H'.ord]=sub_matrix['H'.ord]['A'.ord]=1.0
  sub_matrix['C'.ord]['H'.ord]=sub_matrix['H'.ord]['C'.ord]=1.0
  sub_matrix['T'.ord]['H'.ord]=sub_matrix['H'.ord]['T'.ord]=1.0
  sub_matrix['A'.ord]['V'.ord]=sub_matrix['V'.ord]['A'.ord]=1.0
  sub_matrix['C'.ord]['V'.ord]=sub_matrix['V'.ord]['C'.ord]=1.0
  sub_matrix['G'.ord]['V'.ord]=sub_matrix['V'.ord]['G'.ord]=1.0
  #other
  sub_matrix['$'.ord]['$'.ord]=50.0
  sub_matrix['T'.ord]['U'.ord] = sub_matrix['U'.ord]['T'.ord] = 1.0
  sub_matrix['N'.ord]['N'.ord] = 0.0
  sub_matrix['X'.ord]['-'.ord]=sub_matrix['-'.ord]['X'.ord]=3.0
  ['A','T','G','C'].each do |ch|
    sub_matrix[ch.ord]['*'.ord]=sub_matrix['*'.ord][ch.ord]=1.0
    sub_matrix[ch.ord]['&'.ord]=sub_matrix['&'.ord][ch.ord]=0.7
    sub_matrix[ch.ord]['$'.ord]=sub_matrix['$'.ord][ch.ord]=0.0
    sub_matrix[ch.ord]['.'.ord]=sub_matrix['.'.ord][ch.ord]=-20.0
    sub_matrix[ch.ord]['N'.ord]=sub_matrix['N'.ord][ch.ord]=-3.0
  end
  sub_matrix.each {|column| column.freeze}
  sub_matrix.freeze

  NUCLEOTIDE_MATRIX = sub_matrix

  def self.score_alignment(standard, query)
    sc = 0.0
    0.upto(standard.size - 1) do |i|
      sc += NUCLEOTIDE_MATRIX[standard[i,1].upcase.ord][query[i,1].upcase.ord]
    end
    return sc
  end

  # For a given sequence, find all gaps (continuous sections of gap_char).
  #
  # @param seq [String] The sequence, with gaps marked by gap_char.
  # @param gap_char [Char] The char used to mark gaps in seq. Defaults to '-'.
  # @return [Array<Array<Integer>>] List of gaps (a list of positions).
  def self.make_gap_list(seq, gap_char='-')
    list = []
    cur_ins = nil
    prev_i = nil
    0.upto(seq.size - 1) do |i|
      # 
      if seq[i,1] == gap_char
        # Position is part of gap
        if prev_i and i == prev_i + 1
          # Add current position to existing gap
          cur_ins << i
        else
          # Add current group to list of gaps
          list << cur_ins unless cur_ins.nil? or cur_ins.empty?
          # Start new gap with current position
          cur_ins = [i]
        end
        prev_i = i
      end
    end
    list << cur_ins unless cur_ins.nil? or cur_ins.empty?
    return list
  end

  def self.trim_leading_dashes(standard, query)
    leading_dashes_match = /^(-+)[^-]/.match(standard)
    return if leading_dashes_match.nil?

    leading_dashes = leading_dashes_match[1]
    standard[0, leading_dashes.size] = ''
    query[0, leading_dashes.size] = ''
  end

  def self.trim_trailing_dashes(standard, query)
    trailing_dashes_match = /[^-](-+)$/.match(standard)
    return if trailing_dashes_match.nil?

    trailing_dashes = trailing_dashes_match[1]
    end_of_standard = standard.size - trailing_dashes.size
    standard[end_of_standard, trailing_dashes.size] = ''
    query[end_of_standard, trailing_dashes.size] = ''
  end

  def self.fix_incomplete_edge_codon(query, side=:leading)
    edge_idx = 0
    dash_regex = /^(-+)[^-]/
    incr = 1
    if side != :leading  # fix the trailing edge
      edge_idx = -1
      dash_regex = /[^-](-+)$/
      incr = -1
    end

    if query[edge_idx] == '-'
      dashes = dash_regex.match(query)[1]  # we know there will be a match

      # If the length of the dashes aren't a multiple of 3, turn some
      # of the query characters into dashes to force it to be a full
      # codon of dashes.
      if dashes.size % 3 >= 1
        first_non_dash_idx = dashes.size
        if side != :leading
          first_non_dash_idx = query.size - dashes.size - 1
        end
        query[first_non_dash_idx] = '-'
        if dashes.size % 3 == 1
          query[first_non_dash_idx + incr] = '-'
        end
      end
    end
  end

  def self.merge_insertions_and_deletions_to_fix_of_sequences(
    standard,
    query
  )
    # Merge deletions and insertions until the sequences have a cogent length
    # (i.e. have length divisible by 3).  This helps fix poor insertions near 
    # the start of the sequence.
    raise 'Standard and query should be the same length' if standard.size != query.size
    unless standard.size % 3 == 0
      dex = 0
      while(dex = standard.index(/-/, dex))
        [-1, 1, -2, 2].each do |offset|  # look one base away, then two bases away
          if (dex + offset >= 0) and query[dex + offset] == '-'
            standard[dex] = ''
            query[dex + offset] = ''
            dex = 0
            break
          end
        end
        
        # Stop if the sequences are now a cogent length.
        break if standard.size % 3 == 0
        dex += 1
      end
    end
  end

  # Merges neighbouring gaps to try and force the lengths of all gaps to be a
  # multiple of 3. A list of thresholds defines the maximum distance between
  # the first and the last gap in the list, depending on the number of gaps in
  # the cluster. 
  #
  # @param gaps [Array<Array<Integer>>] List of gaps (a list of positions).
  # @param raise_errors [Boolean] Raise error if clustering fails. 
  #                     Defaults to false.
  # @param thresholds [Array<Integer>] List of distance thresholds in the form:
  #                   [threshold_for_group_of_2, threshold_for_group_of_3, ...]. 
  #                   Defaults to [9, 12].
  # @param trim_distance [Boolean] Disregard inner gap positions in distance.
  def self.cluster_gaps(
    gaps, 
    raise_errors=false, 
    thresholds=[9, 12], 
    trim_distance=false
  )
    # Merge adjacent gaps if they are not a codon-sized gap.
    new_gap_list = []
    gaps.each_with_index do |gap, i|
      # Skip and discard already merged gaps
      next if gap.size == 0
      # Skip and add already valid gaps 
      if gap.size % 3 == 0
        new_gap_list << gap
        next
      end
      # Merge with neighbouring gaps if possible
      did_cluster = false
      thresholds.each.with_index(1) do |threshold, num_other_gaps|
        selection = gaps[i..i+num_other_gaps]
        next unless gaps.size > i + num_other_gaps
        if should_cluster?(selection, threshold, trim_distance)
          new_gap_list << cluster(selection)
          gaps[i..i+num_other_gaps] = [[] * (num_other_gaps + 1)]
          did_cluster = true
          break
        end
      end
      # We can't merge the gaps; either raise an error or meekly proceed.
      unless did_cluster
        if raise_errors
          raise GapMergeError
        else
          # FIXME this behaviour differs between insertions and deletions
          new_gap_list << gap  
        end
      end
    end
    return new_gap_list
  end

  def self.align_gaps_to_frame(gaps, common_gap_locations=nil)
    # Align gaps to codon boundaries, giving preference to common 
    # gap locations if specified.
    # Gaps must be listed in ascending order, i.e. from left to right.
    
    offset = 0  # offset created by previous gaps.
    gaps.each do |gap|
      # See if this gap is close to a common gap location (within 3 amino acids).
      unless common_gap_locations.nil?
        closest_common = common_gap_locations.min do |a, b|
          (3 * a - (gap.first - offset)).abs <=> (3 * b - (gap.first - offset)).abs
        end
        unless closest_common.nil? or (3 * closest_common - (gap.first - offset)).abs > 9
          # Align the gap to this position.
          new_gap = []
          0.upto(gap.size - 1) do |i|
            new_gap << 3 * closest_common + i + offset
          end
          gap.replace(new_gap)
        end
      end
      
      # Align the gap to the nearest appropriate frame.
      # Original comment from Conan: scoring would be good here
      case gap.first % 3
      when 1  # Move back one base
        gap.replace(gap.map(&:pred))
      when 2  # Move forward one base
        gap.replace(gap.map(&:next))
      end
      
      offset += gap.size
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
        if i > seq.size
          seq.insert(-1, '-')
        else
          seq.insert(i, '-')
        end
      end
    end
    return seq
  end

  # Aligns a query sequence to a standard sequence. Assumes standard in the first base.
  #
  # @param standard [String] Standard sequence to align query to.
  # @param query [String] Sequence to align to standard.
  # @param gap_init [Integer] Gap initialization penalty. Defaults to 3.
  # @param gap_extend [Integer] Gap extension penalty. Defaults to 1.
  # @param common_insert_locations [Array<Integer>] List of common insert positions 
  #                                (amino acid positions; starting at base 0). Defaults to [].
  # @param trim [Boolean] Trim trailing and leading insertions/deletions. Defaults to false.
  # @param raise_errors [Boolean] Raise errors when frame alignment fails. Defaults to false.
  # @param prealigned [Boolean] Assume standard and query are already aligned and just run
  #                   corrections and QC. Defaults to false.
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
    if common_insert_locations.nil?
      common_insert_locations = []
    end
    unless prealigned
      standard, query = align_it(standard, query, gap_init, gap_penalty)
    end
    raise "Standard and query should be the same length" unless standard.size == query.size
    
    # Trim leading and trailing dashes if desired.
    if trim
      trim_leading_dashes(standard, query)
      trim_trailing_dashes(standard, query)
      fix_incomplete_edge_codon(query, :leading)
      fix_incomplete_edge_codon(query, :trailing)
    end
    
    merge_insertions_and_deletions_to_fix_of_sequences(standard, query)

    if standard.count('-') % 3 != 0 and raise_errors
      raise "Cannot frame align, #{standard.count('-')} inserted bases not divisible by 3"
    end
    if query.count('-') % 3 != 0 and raise_errors
      raise "Cannot frame align, #{query.count('-')} deleted bases not divisible by 3"
    end
    
    # Build the insert/delete lists.  These lists look like 
    # [[3,4,5], [9], [11,12]]
    insert_list = make_gap_list(standard)
    delete_list = make_gap_list(query)
    
    # Process the insertions.
    unless insert_list.empty?
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
    unless delete_list.empty?
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
    return remove_insertions_from_query(elem[0], elem[1])
  end

  def self.remove_insertions_from_query(standard, query)
    seq = '' + query
    inserts = []
    
    insert_list = []
    0.upto(standard.size - 1) do |i|
      insert_list << i if standard[i,1] == '-'
    end
    
    big_insert_list = []
    if standard.include?('-') # Inserts first
      # First step should be to cluster inserts
      cur_ins = nil
      prev_i = nil
      insert_list.each do |i|
        if prev_i and i == prev_i + 1
          cur_ins << i
        else
          big_insert_list << cur_ins unless cur_ins.nil? or cur_ins.empty?
          cur_ins = [i]
        end
        prev_i = i
      end
      big_insert_list << cur_ins unless cur_ins.nil? or cur_ins.empty?
    end
    
    offset = 0
    big_insert_list.each do |ins|
      ins_seq = ''
      ins.each do |i|
        ins_seq += query[i,1]
      end
      inserts << [((ins.first - offset) / 3), ins_seq]
      offset += ins.size
      ins.each do |i|
        seq[i,1] = '.'
      end
    end

    return [seq.gsub('.',''), inserts]
  end

  private

  # Checks whether the combined length of the gaps is a multiple of three and 
  # whether the distance between the first and the last gap is within the threshold.
  #
  # @param gaps [Array<Array<Integer>>] List of gaps, each represented by a list of positions.
  # @param threshold [Integer] Maximal distance between the first and the last gap.
  # @param trim_distance [Boolean] Disregard inner gap positions in distance calculation.
  # @return [Boolean] The gaps should cluster.
  def self.should_cluster?(gaps, threshold, trim_distance = false)
    return false if gaps.any?(&:nil?)
    cluster_size = gaps.flatten.size
    distance = gaps.last.first - gaps.first.last
    # Remove size of inner gaps from distance
    distance -= gaps[1..-1].flatten.size if trim_distance
    return cluster_size % 3 == 0 && distance < threshold
  end

  # Find the index of the largest gap in the list.
  # 
  # @param gaps [Array<Array<Integer>>] List of gaps (a list of positions).
  # @return [Integer] The index of the largest gap in the list.
  def self.max_gap_index(gaps)
    return (0...gaps.size).max{ |a, b| gaps[a].size <=> gaps[b].size}
  end

  # Merges a list of gaps into a cluster.
  # 
  # @param gaps [Array<Array<Integer>>] List of gaps (a list of positions).
  # @return [Array<Integer>] Merged gap as a list of positions.
  def self.cluster(gaps)
    result = []
    # Find the center gap of the cluster
    if gaps.size % 2 == 1 # Use middle gap as center 
      center = (gaps.size - 1) / 2
    else # Choose the larger of the two middle gaps as center
      center = gaps.size / 2 - 1 + max_gap_index(gaps[(gaps.size / 2 - 1), 2])
    end
    # Compress gaps to the left og the center gap
    if center > 0
      pre = [0, center - 1].max
      # Extend center gap by combined size of all preceding gaps
      result += ((gaps[center].first - gaps[0..pre].flatten.size) .. (gaps[center].first - 1)).to_a
    end
    # Add center gap
    result += gaps[center]
    # Compress gaps to the right of the center gap
    if center < gaps.size - 1
      suc = [center + 1, gaps.size - 1].min
      # Extend center gap by combined size of all subsequent gaps
      result += ((gaps[center].last + 1) .. (gaps[center].last + gaps[suc..gaps.size - 1].flatten.size)).to_a
    end
    return result
  end
end
