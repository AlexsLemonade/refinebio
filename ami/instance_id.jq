# From all of the instances
.Reservations[].Instances[0]
# Take those that have a "Name" tag set to "AMI Template Instance"
| select(
    .Tags[]
    | select(.Key == "Name")
    | select(.Value == "AMI Template Instance")
)
# Then look for the one that is either running or stopped
| select(.State.Name == "running" or .State.Name == "stopped")
# And get its instance id
| .InstanceId
