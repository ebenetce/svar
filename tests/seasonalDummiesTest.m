classdef seasonalDummiesTest < matlab.unittest.TestCase
    %seasonalDummiesTest Tests seasonal dummy construction.

    methods (Test)
        function dropsLastCategoryByDefault(testCase)
            actual = svar.seasonalDummies(14, 12);
            expected = [eye(11); zeros(1, 11); eye(2, 11)];

            testCase.verifyEqual(actual, expected, AbsTol=0);
        end

        function canKeepAllCategories(testCase)
            actual = svar.seasonalDummies(5, 4, DropLast=false);
            expected = [eye(4); 1 0 0 0];

            testCase.verifyEqual(actual, expected, AbsTol=0);
        end

        function dropsSelectedReferenceCategory(testCase)
            actual = svar.seasonalDummies(6, 3, Reference=2);
            expected = [1 0; 0 0; 0 1; 1 0; 0 0; 0 1];

            testCase.verifyEqual(actual, expected, AbsTol=0);
        end

        function rejectsReferenceBeyondPeriod(testCase)
            testCase.verifyError(@() svar.seasonalDummies(12, 4, ...
                Reference=5), "seasonalDummies:InvalidReference");
        end
    end
end
