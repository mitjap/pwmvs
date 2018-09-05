#include "progress.hpp"

#include "exception"

AbstractProgress::AbstractProgress(int total_steps)
    : completed_steps(0)
    , total_steps(0)
{
    configure(total_steps);
}

AbstractProgress::~AbstractProgress()
{

}

void AbstractProgress::increment(int steps)
{
    if (total_steps < 0)
        throw std::logic_error("Invalid step incremet");

    completed_steps = std::min(completed_steps + steps, total_steps);
}

void AbstractProgress::configure(int steps)
{
    if (total_steps < 0)
        throw std::logic_error("Invalid total steps'");

    completed_steps = 0;
    total_steps = steps;
}

int AbstractProgress::completedSteps() const
{
    return completed_steps;
}

int AbstractProgress::totalSteps() const
{
    return total_steps;
}

FloatT AbstractProgress::progress() const
{
    return static_cast<FloatT>(completed_steps) / static_cast<FloatT>(total_steps);
}

ProgressIncrementor::ProgressIncrementor(AbstractProgress *progress, int steps)
    : progress(progress)
    , steps(steps)
{
    // nothing to do
}

ProgressIncrementor::~ProgressIncrementor()
{
    if (progress) progress->increment(steps);
}

ConsoleProgress::ConsoleProgress(int total_steps)
    : AbstractProgress(total_steps)
{
    // nothing to do
}

void ConsoleProgress::increment(int steps)
{
    AbstractProgress::increment(steps);

    std::cout << "Progress: " << completedSteps() << "/" << totalSteps() << "\n" << std::endl;
}

void ConsoleProgress::configure(int total_steps)
{
    int old_total_steps = totalSteps();
    AbstractProgress::configure(total_steps);

    if (old_total_steps != totalSteps()) {
        increment(0); // Just to print out the progress
    }
}
