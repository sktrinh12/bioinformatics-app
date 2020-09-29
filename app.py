import streamlit as st
# from Bio.Seq import Seq
from Bio import SeqIO
from functions import *
import neatbio.sequtils as utils
from collections import Counter
import io


# disable warning of deprecation
st.set_option('deprecation.showfileUploaderEncoding', False)


def main():
    """A simple bioinformatics app"""
    st.title("Simple Bioinformatics App")
    menu = ['Intro', 'DNA Sequence', 'DotPlot', 'About']

    choice = st.sidebar.selectbox('Select Activity', menu)
    if choice == 'Intro':
        st.subheader('Intro to BioInformatics')
    elif choice == 'DNA Sequence':
        st.subheader('DNA Sequence Analysis')
        seq_file = st.file_uploader('Upload FASTA file', type=['fasta', 'fa',
                                                               'txt', 'fna'])

        if seq_file is not None:
            # to no longer autodetect the file's encoding.
            # This means that all files will be returned as binary buffers.
            # thus must wrap in textio wrapper
            # seq_file_io = io.TextIOWrapper(seq_file)
            dna_record = SeqIO.read(seq_file, 'fasta')
            dna_seq = dna_record.seq

            details = st.radio('Details', ('Description','Sequence'))
            if details == 'Description':
                st.write(dna_record.description)
            elif details == 'Sequence':
                st.write(dna_seq)


            # nucleotide frequencies
            st.subheader('Nucleotide Frequency')
            dna_freq = Counter(dna_seq)
            st.write(dna_freq)
            adenine_colour = st.beta_color_picker('Adenine Colour')
            guanine_colour = st.beta_color_picker('Guanine Colour')
            cytosine_colour = st.beta_color_picker('Cytosine Colour')
            thymine_colour = st.beta_color_picker('Thymine Colour')

            if st.button('Plot Frequency'):
                fig, ax = plt.subplots()
                barlist = ax.bar(dna_freq.keys(), dna_freq.values())
                barlist[3].set_color(cytosine_colour)
                barlist[0].set_color(guanine_colour)
                barlist[2].set_color(adenine_colour)
                barlist[1].set_color(thymine_colour)
                st.pyplot(fig)

            st.subheader('DNA Composition')
            gc_score = utils.gc_content(str(dna_seq))
            at_score = utils.at_content(str(dna_seq))
            st.write({'GC Content' : gc_score, 'AT Content' : at_score})
            st.json({'GC Content' : gc_score, 'AT Content' : at_score})

            # nucleotide count
            nt_count = st.text_input('Enter Nucleotide Here', 'Type Nucleotide Alphabet')
            st.write('Number of {} Nucleotide is ::{}'.format(nt_count,\
                                                              str(dna_seq).count(nt_count)))

            # protein synthesis
            st.subheader('Protein Synthesis')
            p1 = dna_seq.translate()
            aa_freq = Counter(str(p1))

            if st.checkbox('Transcription'):
                st.write(dna_seq.transcribe())

            elif st.checkbox('Translate'):
                st.write(p1)

            elif st.checkbox('Complement'):
                st.write(dna_seq.complement())

            elif st.checkbox('AA Frequency'):
                st.write(aa_freq)

            # top most common amino acid
            elif st.checkbox('Plot AA Frequency'):
                fig, ax = plt.subplots()
                aa_colour = st.beta_color_picker('Pick an amino acid colour')
                barlist = ax.bar(aa_freq.keys(), aa_freq.values(),
                                 color=aa_colour)
                st.pyplot(fig)

            elif st.checkbox('Full amino acid name'):
                aa_name = str(p1).replace('*', '')
                aa3 = utils.convert_1to3(aa_name)
                st.write(aa_name)
                st.write('='*30)
                st.write(aa3)

                st.write('='*30)
                st.write(utils.get_acid_name(aa3))




    elif choice == 'DotPlot':
        st.subheader('Generate Dot Plot for two Sequences')

        seq_file1 = st.file_uploader('Upload 1st FASTA file', type=['fasta', 'fa',
                                                               'txt', 'fna'])

        seq_file2 = st.file_uploader('Upload 2nd FASTA file', type=['fasta', 'fa',
                                                               'txt', 'fna'])
        if seq_file1 and seq_file2 is not None:
            dna_record1 = SeqIO.read(seq_file1, 'fasta')
            dna_record2 = SeqIO.read(seq_file2, 'fasta')
            dna_seq1 = dna_record1.seq
            dna_seq2 = dna_record2.seq

            details = st.radio('Details', ('Description','Sequence'))
            if details == 'Description':
                st.write(dna_record1.description)
                st.write('='*50)
                st.write(dna_record2.description)
            elif details == 'Sequence':
                st.write(dna_record1.seq)
                st.write('='*50)
                st.write(dna_record2.seq)

            # (label, min_value, max_value, value)
            custom_limit = st.number_input('Select max number of nucleotides',
                                           10, 400, 50)
            if st.button('Dot Plot'):
                seq1 = join_str(dna_seq1, custom_limit)
                seq2 = join_str(dna_seq2, custom_limit)
                numeric_arr1 = assign_numeric_array(seq1)
                numeric_arr2 = assign_numeric_array(seq2)
                compr_array = generate_comparison_array(numeric_arr1,
                                                        numeric_arr2)
                fig, ax = plot_dotplot(compr_array, seq1, seq2)
                st.write('Comparing the first {} nucleotides of the two sequences'.format(custom_limit))
                st.pyplot(fig)


    elif choice == 'About':
        st.subheader('About')




if __name__ == '__main__':
    main()
