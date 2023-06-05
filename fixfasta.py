import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='addes > to the fasta') 
    parser.add_argument('file', type=str, help='fixes file')
    args = parser.parse_args()

    with open(args.file, 'r') as reader:
        with open(args.file + '.fixed', 'w') as writer:
            line = reader.readline()
            while line != '':
                if '@' in line:
                    writer.write('>' + line)
                else:
                    writer.write(line)
                line = reader.readline()

print('done')
            