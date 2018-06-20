# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/heuristics/Band.cpp \
../src/heuristics/BandCalculator.cpp \
../src/heuristics/Node.cpp \

OBJS += \
./src/heuristics/Band.o \
./src/heuristics/BandCalculator.o \
./src/heuristics/Node.o \

CPP_DEPS += \
./src/heuristics/Band.d \
./src/heuristics/BandCalculator.d \
./src/heuristics/Node.d \


# Each subdirectory must supply rules for building sources it contributes
src/heuristics/%.o: ../src/heuristics/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	$(CC) $(INC_PAR) $(CPPFLAGS) -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


