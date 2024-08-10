from datetime import datetime, timedelta, date


def THP_BIOID_Time(PMA_Treat_Time, Day=str(date.today())):
    """
    This function calculates the time in which steps should occur for the THP-1 BioID proximity
    localizaiton assay experiment.

    Args:
        PMA_Treat_Time (str): The treatment start time given in Military time 24-hour format (HH:MM)
        Day (str): The current day of the treatment (YYY-MM-DD)

    Returns:
        PMA_Media_Exchange: The time at which you will need to change out the PMA treated RPMI Media for DMEM
        Dox_induction: The time you will need to start the Dox induction
        Biotin_Treatment: The time you will need to treat the cells with biotin
        Lysis: The time in which you will need to start lysing the cells
    """

    Day_of_week = {
        0: 'Monday',
        1: 'Tuesday',
        2: 'Wednesday', 
        3: 'Thursday', 
        4: 'Friday',
        5: 'Saturday', 
        6: 'Sunday'
    }

    PMA_Treat_Time = datetime.strptime(PMA_Treat_Time, "%H:%M")
    Day = datetime.strptime(Day, "%Y-%m-%d")
    Exp_Start = datetime.combine(Day, PMA_Treat_Time.time())
    PMA_Media_Exchange = Exp_Start + timedelta(hours=24)
    Dox_Induction = PMA_Media_Exchange + timedelta(days=2)
    Biotin_Treatment = Dox_Induction + timedelta(hours=32)
    Lysis = Dox_Induction + timedelta(days=2)

    print(f"\033[33mYou will PMA treat your THP-1 Cells on: \033[0m {Day_of_week[Exp_Start.weekday()]}, {Exp_Start}")
    print(f"\033[33mPMA_Media_Exchange: \033[0m {Day_of_week[PMA_Media_Exchange.weekday()]}, {PMA_Media_Exchange}")
    print(f"\033[33mDox_Induction: \033[0m {Day_of_week[Dox_Induction.weekday()]}, {Dox_Induction}")
    print(f"\033[33mBiotin_Treatment: \033[0m {Day_of_week[Biotin_Treatment.weekday()]}, {Biotin_Treatment}")
    print(f"\033[33mLysis: \033[0m {Day_of_week[Lysis.weekday()]}, {Lysis}")

    return PMA_Media_Exchange, Dox_Induction, Biotin_Treatment, Lysis

# Determine whether I have already started the experiment
experiment_input = input('\033[35mWill you start your experiment today or another day?\033[0m')

# Determines whether I will use todays date automatically or another date
if experiment_input.lower() == 'today':
    THP_BIOID_Time(input('\033[35mWhat time will you start your experiment (HH:MM): \033[0m'))
else:
    THP_BIOID_Time(input('\033[35mWhat time will you start your experiment (HH:MM): \033[0m'), \
                   input('\033[35mWhat day will you start your experiment (YYYY-MM-DD): \033[0m'))


    