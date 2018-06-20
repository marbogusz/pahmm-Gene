# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/pHMMg.cpp 

OBJS += \
./src/pHMMg.o 

CPP_DEPS += \
./src/pHMMg.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	$(CC) $(INC_PAR) $(CPPFLAGS) -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


