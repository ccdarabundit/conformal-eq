cmake_minimum_required(VERSION 3.15)
project(ConformalEQ VERSION 0.9)

set(CMAKE_CXX_STANDARD 17)
add_subdirectory(modules)

juce_add_plugin(ConformalEQ
	COMPANY_NAME Champion Sound Systems
	PLUGIN_MANUFACTURER_CODE CSS
	PLUGIN_CODE cfeq
	FORMATS AU VST3 Standalone
	ProductName "Conformal EQ")

juce_generate_juce_header(ConformalEQ)

add_subdirectory(src)
include_subdirectory(src)

target_compile_definitions(ConformalEQ PUBLIC
    JUCE_VST3_CAN_REPLACE_VST2=0
)

target_link_libraries(ConformalEQ PRIVATE 
		juce_plugin_modules)