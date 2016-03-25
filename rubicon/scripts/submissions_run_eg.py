# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

from argparse import ArgumentParser

from rubicon.processors.process_submissions_eg import SubmissionProcessorEG


def go_eg_submissions():
    m_description = 'This program is used to pull jobs from the Submissions ' \
                    'database, create FireWorks workflows from those ' \
                    'submissions, and then monitor all previous submissions ' \
                    'for updates to state (so that the submission ' \
                    'database can be updated)'

    parser = ArgumentParser(description=m_description)
    parser.add_argument('--sleep', help='sleep time between loops',
                        default=None, type=int)
    parser.add_argument('--infinite', help='loop infinite times',
                        action='store_true')
    args = parser.parse_args()

    sp = SubmissionProcessorEG.auto_load()
    sp.run(args.sleep, args.infinite)


if __name__ == '__main__':
    go_eg_submissions()
