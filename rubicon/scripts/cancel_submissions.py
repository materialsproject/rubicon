import argparse
import datetime
from rubicon.submission.submission_mongo_eg import SubmissionMongoAdapterEG

__author__ = 'xiaohuiqu'



def main():
    parser = argparse.ArgumentParser(
        description="Cancel the submitted jobs within a specific duration, default to 6 hours")
    parser.add_argument("-d", "--days", dest="days", type=int, default=0,
                        help="the duration of days")
    parser.add_argument("-r", "--hours", dest="hours", type=int, default=6,
                        help="the duration of hours")
    parser.add_argument("-m", "--minutes", dest="minutes", type=int, default=0,
                        help="the duration of minutes")
    options = parser.parse_args()
    t_delta = datetime.timedelta(days=options.days, hours=options.hours, minutes=options.minutes)
    t_now = datetime.datetime.utcnow()
    t_submission = t_now - t_delta
    sma = SubmissionMongoAdapterEG.auto_load()
    ids = sma.get_submission_ids_after(t_submission)
    for i in ids:
        sma.set_job_state_to_cancel(i)

if __name__ == '__main__':
    main()