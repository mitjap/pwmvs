#pragma once

#include "types.hpp"

class AbstractProgress
{
public:
    //AbstractProgress();
    AbstractProgress(int total_steps = 0);
    virtual ~AbstractProgress();

    virtual void increment(int steps = 1);
    virtual void configure(int total_steps);

    int completedSteps() const;
    int totalSteps() const;
    FloatT progress() const;

private:
    int total_steps;
    int completed_steps;
};

class ProgressIncrementor
{
public:
    ProgressIncrementor(AbstractProgress *progress, int steps = 1);
    virtual ~ProgressIncrementor();

private:
    AbstractProgress *progress;
    int steps;
};

class ConsoleProgress : public AbstractProgress
{
public:
    ConsoleProgress(int total_steps = 0);
    virtual void increment(int steps = 1);
    virtual void configure(int total_steps);
};
